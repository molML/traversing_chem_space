import numpy as np
import torch
from torch import Tensor, tensor
from torch_geometric.loader import DataLoader as pyg_DataLoader
from torch.utils.data import TensorDataset, DataLoader
from tqdm.auto import trange, tqdm
import pyro
from pyro.infer.autoguide import AutoDiagonalNormal
from pyro.infer import SVI, Trace_ELBO, Predictive

from active_learning.nn.GNNs import GCN
from active_learning.nn.bayesian import BNN, BayesianGCNModel
from active_learning.hyperopt import optimize_hyperparameters


class GCNEnsemble:
    """ Ensemble of GCNs"""
    def __init__(self, ensemble_size: int = 10, seed: int = 0, **kwargs) -> None:
        self.ensemble_size = ensemble_size
        self.seed = seed
        rng = np.random.default_rng(seed=seed)
        self.seeds = rng.integers(0, 1000, ensemble_size)
        self.models = {i: GCN(seed=s, **kwargs) for i, s in enumerate(self.seeds)}

    def optimize_hyperparameters(self, x, y, **kwargs):
        best_hypers = optimize_hyperparameters(x, y, architecture='gcn', **kwargs)
        # re-init model wrapper with optimal hyperparameters
        self.__init__(ensemble_size=self.ensemble_size, seed=self.seed, **best_hypers)

    def train(self, x, y, epochs: int = None, verbose: bool = False, **kwargs) -> None:
        for i, m in self.models.items():
            m.train(x, y, epochs, verbose=verbose, **kwargs)

    def predict(self, x, **kwargs) -> (Tensor, Tensor, Tensor):

        y_hat = torch.stack([m.predict(x, **kwargs) for m in self.models.values()]).T
        # to logits
        y_hat_mu = torch.mean(y_hat.float(), 1)
        # how close to the decision boundary is the predicted value. Closer is more uncertain, with 0.5 being the max
        y_hat_sigma = abs((y_hat_mu > 0.5) * 1 - y_hat_mu)

        return y_hat, y_hat_mu, y_hat_sigma

    def __getitem__(self, item):
        return self.models[item]

    def __repr__(self) -> str:
        return f"Ensemble of {self.ensemble_size} GCN Classifiers"


class BayesianNN:
    """ Bayesian Neural Network wrapper with train() and predict() functions. """

    def __init__(self, seed: int = 42, to_gpu: bool = False, **kwargs):

        # Define some vars and seed random state
        self.device = torch.device("cuda:0" if torch.cuda.is_available() and to_gpu else "cpu")
        self.train_losses = []
        self.to_gpu = to_gpu
        self.epoch = 0
        self.seed = seed
        pyro.set_rng_seed(seed)

        # Init model
        self.model = BNN(to_gpu=to_gpu, **kwargs).to(self.device)
        self.lr = self.model.lr
        self.epochs = self.model.epochs

        # Init Guide model
        self.guide = AutoDiagonalNormal(self.model)
        self.guide = self.guide.to(self.device)
        # Init optimizer
        adam = pyro.optim.Adam({"lr": self.lr})
        # Stochastic Variational Inference
        self.svi = SVI(self.model, self.guide, adam, loss=Trace_ELBO())

    def train(self, x: np.ndarray, y: np.ndarray, epochs: int = None, batch_size: int = 256) -> None:

        graph = True if type(x) is not np.ndarray else False

        # Convert numpy to Torch
        if graph:
            data_loader = pyg_DataLoader(x, batch_size=batch_size)
        else:
            data_loader = DataLoader(TensorDataset(Tensor(x), Tensor(y)), batch_size=batch_size)

        # Training loop
        pyro.clear_param_store()
        bar = trange(self.epochs if epochs is None else epochs)
        for epoch in bar:

            running_loss, n_samples = 0.0, 0
            for batch in data_loader:
                if graph:
                    batch = batch.to(self.device)
                    # ELBO gradient and add loss to running loss
                    running_loss += self.svi.step(x=batch.x.float(), edge_index=batch.edge_index,
                                                  batch=batch.batch, y=batch.y)
                    n_samples += len(batch)
                else:
                    x, y = batch[0].to(self.device), batch[1].to(self.device)
                    # ELBO gradient and add loss to running loss
                    running_loss += self.svi.step(x, y)
                    n_samples += x.shape[0]

            loss = running_loss / n_samples
            self.train_losses.append(loss)
            bar.set_postfix(loss=f'{loss:.4f}')

    def predict(self, x, num_samples: int = 500, batch_size: int = 256) -> (tensor, tensor, tensor):

        # Construct predictive distribution
        predictive = Predictive(self.model, guide=self.guide, num_samples=num_samples, return_sites=("obs", "_RETURN"))
        y_hat = []

        # Convert numpy to Torch data loader and predict
        graph = True if type(x) is not np.ndarray else False
        if graph:
            data_loader = pyg_DataLoader(x, batch_size=batch_size)
        else:
            data_loader = DataLoader(TensorDataset(Tensor(x)), batch_size=batch_size)

        for batch in tqdm(data_loader, 'Sampling predictive distribution'):
            if graph:
                batch = batch.to(self.device)
                # Reshape if needed
                samples = predictive(x=batch.x.float(), edge_index=batch.edge_index, batch=batch.batch)
                y_hat.append(samples['obs'].T)

            else:
                x = batch[0].to(self.device)
                # Reshape if needed
                samples = predictive(x.unsqueeze(0) if len(x.size()) == 1 else x)
                y_hat.append(samples['obs'].T)

        y_hat = torch.cat(y_hat, dim=0)
        # to logits
        y_hat_mu = torch.mean(y_hat.float(), 1)
        # how close to the decision boundary is the predicted value. Closer is more uncertain, with 0.5 being the max
        y_hat_sigma = abs((y_hat_mu > 0.5) * 1 - y_hat_mu)

        return y_hat, y_hat_mu, y_hat_sigma

    def optimize_hyperparameters(self, x, y, **kwargs):
        best_hypers = optimize_hyperparameters(x, y, architecture='bnn', **kwargs)
        # re-init model wrapper with optimal hypers
        self.__init__(seed=self.seed, to_gpu=self.to_gpu, **best_hypers)


class BayesianGCN(BayesianNN):
    def __init__(self, seed: int = 42, to_gpu: bool = False, lr=1e-3, epochs: int = 10000, **kwargs):
        # Define some vars and seed random state
        self.lr = lr
        self.train_losses = []
        self.epochs = epochs
        self.epoch = 0
        self.seed = seed
        self.device = torch.device("cuda:0" if torch.cuda.is_available() and to_gpu else "cpu")
        pyro.set_rng_seed(seed)

        # Init model
        self.model = BayesianGCNModel(to_gpu=to_gpu, **kwargs).to(self.device)

        # Init Guide model
        self.guide = AutoDiagonalNormal(self.model)
        self.guide = self.guide.to(self.device)
        # Init optimizer
        adam = pyro.optim.Adam({"lr": lr})
        # Stochastic Variational Inference
        self.svi = SVI(self.model, self.guide, adam, loss=Trace_ELBO())


# from active_learning.data_prep import MasterDataset
# from active_learning.utils import Evaluate
# ds_test = MasterDataset('test', representation='graph')
# x_train, y_train, smiles_train = ds_test[range(10000)]
#
# model = GCNEnsemble(hidden_channels=512)
# model.train(x_train, y_train, verbose=True, epochs=100)
# y_hat, y_hat_mu, y_hat_sigma = model.predict(x_train)
#
# eval = Evaluate()
# eval.eval(y_hat_mu, y_train)
# eval

#
# ds_test = MasterDataset('test', representation='ecfp')
# x_train, y_train, smiles_train = ds_test[range(1000)]
# #
# model = BayesianNN()
# model.train(x_train, y_train, epochs=1000)
# y_hat, y_hat_mu, y_hat_sigma = model.predict(x_train)

