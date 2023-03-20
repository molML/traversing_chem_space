
import numpy as np
import torch
import torch.nn.functional as F
from torch import Tensor
from torch.nn import Linear
from torch_geometric.nn import GCNConv, global_add_pool, GraphNorm
from torch_geometric.loader import DataLoader
from warnings import warn
from tqdm.auto import trange, tqdm
from xgboost import XGBClassifier
import warnings
warnings.filterwarnings("ignore")


class Model:
    def __init__(self, architecture: str = 'xgb', seed: int = 0, **kwargs):
        self.architecture = architecture
        if architecture == 'xgb':
            self.model = XGBoostEnsemble(seed=seed, **kwargs)
        elif architecture == 'gnn':
            self.model = GCNEnsemble(seed=seed, **kwargs)
        elif architecture == 'bnn':
            raise NotImplementedError

    def optimize_hyperparameters(self, **kwargs):
        raise NotImplementedError

    def train(self, x, y, **kwargs):
        self.model.train(x, y, **kwargs)

    def test(self, x, y, **kwargs):
        raise NotImplementedError

    def predict(self, x, **kwargs):
        return self.model.predict(x, **kwargs)

    def __repr__(self):
        return f"{self.architecture} model wrapper"

    def __call__(self, x, **kwargs):
        return self.predict(x, **kwargs)


class XGBoostEnsemble:
    """ Ensemble of n XGBoost regressors, seeded differently """
    def __init__(self, ensemble_size: int = 1, seed: int = 0, **kwargs) -> None:
        self.ensemble_size = ensemble_size
        rng = np.random.default_rng(seed=seed)
        seeds = rng.integers(0, 1000, ensemble_size)
        self.models = {i: XGBClassifier(random_state=s, **kwargs) for i, s in enumerate(seeds)}

    def predict(self, x: np.ndarray, invert_classes: bool = True) -> (np.ndarray, np.ndarray, np.ndarray):
        y_hat_proba = np.array([model.predict_proba(x) for i, model in self.models.items()])

        if invert_classes:
            y_hat_proba = (y_hat_proba * -1) + 1

        y_hat, y_hat_uncertainty = y_hat_proba[:, :, 0], y_hat_proba[:, :, 1]
        mu = np.mean(y_hat, axis=0)
        sigma = np.mean(y_hat_uncertainty, axis=0)

        return y_hat, mu, sigma

    def train(self, x: np.ndarray, y: np.ndarray, **kwargs) -> None:
        for i, m in self.models.items():
            self.models[i] = m.fit(x, y)

    def __getitem__(self, item):
        return self.models[item]

    def __repr__(self) -> str:
        return f"Ensemble of {self.ensemble_size} XGBoost Classifiers"


class GCNEnsemble:
    def __init__(self, ensemble_size: int = 10, seed: int = 0, **kwargs):
        self.ensemble_size = ensemble_size
        rng = np.random.default_rng(seed=seed)

        seeds = rng.integers(0, 1000, ensemble_size)
        self.models = {i: GCN(seed=s, **kwargs) for i, s in enumerate(seeds)}

    def train(self, dataloader: DataLoader, epochs: int = None):
        for i, m in self.models.items():
            m.train(dataloader, epochs)

    def predict(self, dataloader):
        y_hat = torch.stack([m.predict(dataloader) for m in self.models.values()]).T
        mu = torch.mean(y_hat, dim=1)
        sigma = torch.std(y_hat, dim=1)

        return y_hat, mu, sigma

    def __getitem__(self, item):
        return self.models[item]

    def __repr__(self) -> str:
        return f"Ensemble of {self.ensemble_size} GCN Classifiers"


class BayesianGNN:
    def __init__(self):
        raise NotImplementedError


class GCNModel(torch.nn.Module):
    def __init__(self, in_channels: int = 59, hidden_channels: int = 256, num_conv_layers: int = 3, lr: float = 0.001,
                 epochs: int = 100):
        super(GCNModel, self).__init__()
        self.lr = lr
        self.epochs = epochs
        self.in_channels = in_channels
        self.hidden_channels = hidden_channels
        self.num_conv_layers = num_conv_layers

        self.lin_0 = Linear(in_channels, hidden_channels)
        self.lin_1 = Linear(hidden_channels, hidden_channels)
        self.lin_out = Linear(hidden_channels, 1)

        self.convs, self.norms = torch.nn.ModuleList(), torch.nn.ModuleList()
        for _ in range(num_conv_layers):
            self.convs.append(GCNConv(hidden_channels, hidden_channels))
            self.norm = GraphNorm(hidden_channels)

    def reset_parameters(self):
        self.lin_0.reset_parameters()
        self.lin_1.reset_parameters()
        self.lin_out.reset_parameters()
        for conv, norm in zip(self.convs, self.norms):
            conv.reset_parameters()
            norm.reset_parameters()

    def forward(self, x: Tensor, edge_index: Tensor, batch: Tensor) -> Tensor:
        # Atom Embedding:
        x = F.elu(self.lin_0(x))

        # Graph convolutions
        for conv, norm in zip(self.convs, self.norms):
            x = conv(x, edge_index)
            x = norm(x, batch)
            x = F.relu(x)

        # Perform global pooling
        x = global_add_pool(x, batch)

        # Prediction head
        x = F.relu(self.lin_1(x))
        x = self.lin_out(x)

        return torch.sigmoid(x).squeeze()


class GCN:
    def __init__(self, seed: int = 42, **kwargs) -> None:

        torch.manual_seed(seed)
        self.model = GCNModel(**kwargs)
        self.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        self.loss_fn = torch.nn.BCELoss()

        # Move the whole model to the gpu
        self.model = self.model.to(self.device)
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=self.model.lr)

        self.train_loss = []
        self.epochs, self.epoch = self.model.epochs, 0

    def train(self, dataloader: DataLoader, epochs: int = None) -> None:

        bar = trange(self.epochs if epochs is None else epochs)
        for epoch in bar:
            running_loss = 0
            items = 0

            for idx, batch in enumerate(dataloader):
                batch.to(self.device)
                self.optimizer.zero_grad()

                # Forward pass
                y_hat = self.model(batch.x.float(), batch.edge_index, batch.batch)
                loss = self.loss_fn(y_hat, batch.y.float())

                if not loss > 0:
                    warn(f"Failed forward pass on batch {idx}, y_hat = {y_hat}", category=RuntimeWarning)
                else:
                    loss.backward()
                    self.optimizer.step()

                    running_loss += loss.item() * len(batch.y)
                    items += len(batch.y)

            epoch_loss = running_loss / items
            bar.set_postfix(loss=f'{epoch_loss:.4f}')
            self.train_loss.append(epoch_loss)
            self.epoch += 1

    def predict(self, dataloader: DataLoader) -> Tensor:
        """ Predict

        :param dataloader: Torch geometric data loader with data
        :return: A 1D-tensors
        """
        y_hats = torch.tensor([]).to(self.device)
        with torch.no_grad():
            for batch in dataloader:
                batch.to(self.device)
                y_hat = self.model(batch.x.float(), batch.edge_index, batch.batch)
                y_hats = torch.cat((y_hats, y_hat), 0)

        return y_hats


# from active_learning.utils import Evaluate
# from active_learning.data_prep import MasterDataset
# ds_test = MasterDataset('test', representation='graph')
#
# graphs, y, smiles = ds_test[range(10000)]
#
# dataloader = DataLoader(graphs, batch_size=256)
# model = GCN(in_channels=59, epochs=250)
# model.train(dataloader)
# y_hat = model.predict(dataloader)
#
# eval = Evaluate()
# eval.eval(y_hat, y)
# print(eval.to_dataframe())
#
#
# model = GCNEnsemble(ensemble_size=3, epochs=50)
# model.train(dataloader)
# model.models
#
#




