
import numpy as np
import torch
from torch import Tensor, tensor
from torch_geometric.loader import DataLoader as pyg_DataLoader
from torch.utils.data import TensorDataset, DataLoader
from tqdm.auto import trange, tqdm
from torch.nn import functional as F
from warnings import warn
from copy import deepcopy
from torch.nn import Linear
from torch_geometric.nn import GCNConv, global_add_pool, GraphNorm
from warnings import warn
from tqdm.auto import trange
from active_learning.hyperopt import optimize_hyperparameters
import math


class MLP(torch.nn.Module):
    def __init__(self, in_feats: int = 1024, n_hidden: int = 1024, n_out: int = 2, n_layers: int = 3, seed: int = 42,
                 lr: float = 1e-4, epochs: int = 100, anchored: bool = True, l2_lambda: float = 1e-5):
        super().__init__()
        self.seed, self.lr, self.l2_lambda, self.epochs, self.anchored = seed, lr, l2_lambda, epochs, anchored
        torch.manual_seed(seed)

        self.fc = torch.nn.ModuleList()
        for i in range(n_layers):
            self.fc.append(torch.nn.Linear(in_feats if i == 0 else n_hidden, n_hidden))
        self.out = torch.nn.Linear(n_hidden, n_out)

    def reset_parameters(self):
        for lin in self.fc:
            lin.reset_parameters()
        self.out.reset_parameters()

    def forward(self, x: Tensor) -> Tensor:
        for lin in self.fc:
            x = F.relu(lin(x))
        x = self.out(x)
        x = F.log_softmax(x, 1)

        return x


class GCN(torch.nn.Module):
    def __init__(self, in_feats: int = 59, n_hidden: int = 1024, num_conv_layers: int = 3, lr: float = 0.0001,
                 epochs: int = 100, n_out: int = 2, n_layers: int = 3, seed: int = 42, anchored: bool = True,
                 l2_lambda: float = 1e-5):

        super().__init__()
        self.seed, self.lr, self.l2_lambda, self.epochs, self.anchored = seed, lr, l2_lambda, epochs, anchored

        self.atom_embedding = torch.nn.Linear(in_feats, n_hidden)
        self.fc = torch.nn.ModuleList()
        for i in range(n_layers):
            self.fc.append(torch.nn.Linear(n_hidden, n_hidden))
        self.out = torch.nn.Linear(n_hidden, n_out)

        self.convs = torch.nn.ModuleList()
        # self.norms = torch.nn.ModuleList()
        for _ in range(num_conv_layers):
            self.convs.append(GCNConv(n_hidden, n_hidden))
            # self.norm = GraphNorm(n_hidden)

    def reset_parameters(self):
        self.atom_embedding.reset_parameters()
        for conv in self.convs:
            conv.reset_parameters()
            # norm.reset_parameters()
        for lin in self.fc:
            lin.reset_parameters()
        self.out.reset_parameters()

    def forward(self, x: Tensor, edge_index: Tensor, batch: Tensor) -> Tensor:
        # Atom Embedding:
        x = F.relu(self.atom_embedding(x))

        # Graph convolutions
        for conv in self.convs:
            x = conv(x, edge_index)
            x = F.relu(x)

        # Perform global pooling by sum pooling
        x = global_add_pool(x, batch)

        for lin in self.fc:
            x = F.relu(lin(x))
        x = self.out(x)
        x = F.log_softmax(x, 1)

        return x


class Model(torch.nn.Module):
    def __init__(self, architecture: str, class_weights: list = None,  **kwargs):
        super().__init__()
        assert architecture in ['gcn', 'mlp']
        self.architecture = architecture
        self.model = MLP(**kwargs) if architecture == 'mlp' else GCN(**kwargs)

        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        if class_weights is not None:
            class_weights = class_weights.float().to(self.device)
        self.loss_fn = torch.nn.NLLLoss(weight=class_weights)

        # Move the whole model to the gpu
        self.model = self.model.to(self.device)
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=self.model.lr)

        # Save initial weights in the model for the anchored regularization and move them to the gpu
        if self.model.anchored:
            self.model.anchor_weights = deepcopy({i: j for i, j in self.model.named_parameters()})
            self.model.anchor_weights = {i: j.to(self.device) for i, j in self.model.anchor_weights.items()}

        self.train_loss = []
        self.epochs, self.epoch = self.model.epochs, 0

    def train(self, dataloader: DataLoader, epochs: int = None, verbose: bool = True) -> None:

        bar = trange(self.epochs if epochs is None else epochs, disable=not verbose)

        for _ in bar:
            running_loss = 0
            items = 0

            for idx, batch in enumerate(dataloader):

                self.optimizer.zero_grad()

                if self.architecture == 'gcn':
                    batch.to(self.device)
                    y = batch.y
                    y_hat = self.model(batch.x.float(), batch.edge_index, batch.batch)
                else:
                    x, y = batch[0].to(self.device), batch[1].to(self.device)
                    y_hat = self.model(x)

                if len(y_hat) == 0:
                    y_hat = y_hat.unsqueeze(0)
                loss = self.loss_fn(y_hat, y.squeeze())

                if self.model.anchored:
                    # Calculate the total anchored L2 loss
                    l2_loss = 0
                    for param_name, params in self.model.named_parameters():
                        anchored_param = self.model.anchor_weights[param_name]

                        l2_loss += (self.model.l2_lambda / len(y)) * torch.mul(params - anchored_param,
                                                                               params - anchored_param).sum()

                    # Add anchored loss to regular loss according to Pearce et al. (2018)
                    loss = loss + l2_loss

                if not loss >= 0:
                    warn(f"Failed forward pass on batch {idx}, y_hat = {y_hat}", category=RuntimeWarning)
                else:
                    loss.backward()
                    self.optimizer.step()

                    running_loss += loss.item() * len(y)
                    items += len(y)

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
                if self.architecture == 'gcn':
                    batch.to(self.device)
                    y_hat = self.model(batch.x.float(), batch.edge_index, batch.batch)
                else:
                    x = batch[0].to(self.device)
                    y_hat = self.model(x)
                if len(y_hat) == 0:
                    y_hat = y_hat.unsqueeze(0)
                y_hats = torch.cat((y_hats, y_hat), 0)

        return y_hats


class Ensemble(torch.nn.Module):
    """ Ensemble of GCNs"""
    def __init__(self, ensemble_size: int = 10, seed: int = 0, **kwargs) -> None:
        self.ensemble_size = ensemble_size
        self.seed = seed
        rng = np.random.default_rng(seed=seed)
        self.seeds = rng.integers(0, 1000, ensemble_size)
        self.models = {i: Model(seed=s, **kwargs) for i, s in enumerate(self.seeds)}

    def optimize_hyperparameters(self, x, y: DataLoader, **kwargs):
        # raise NotImplementedError
        best_hypers = optimize_hyperparameters(x, y, architecture='gcn', **kwargs)
        # # re-init model wrapper with optimal hyperparameters
        self.__init__(ensemble_size=self.ensemble_size, seed=self.seed, **best_hypers)

    def train(self, dataloader: DataLoader, **kwargs) -> None:
        for i, m in self.models.items():
            m.train(dataloader, **kwargs)

    def predict(self, dataloader, **kwargs) -> Tensor:
        """ logits_N_K_C = [N, num_inference_samples, num_classes] """
        logits_N_K_C = torch.stack([m.predict(dataloader) for m in self.models.values()], 1)

        return logits_N_K_C

    def __getitem__(self, item):
        return self.models[item]

    def __repr__(self) -> str:
        return f"Ensemble of {self.ensemble_size} Classifiers"
