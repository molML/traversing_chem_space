"""

This script contains all models:

    - MLP: a simple feed forward multi-layer perceptron. Supports weight anchoring - Pearce et al. (2018)
    - GCN: a simple graph convolutional NN - Kipf & Welling (2016). Supports weight anchoring - Pearce et al. (2018)
    - Model: A wrapper class that contains a train(), and predict() loop
    - Ensemble: Class that ensembles n Model classes. Contains a train() method and an predict() method that outputs
        logits_N_K_C, defined as [N, num_inference_samples, num_classes]. Also has an optimize_hyperparameters() method.

    Author: Derek van Tilborg, Eindhoven University of Technology, May 2023

"""

import numpy as np
import torch
from torch import Tensor
from torch.utils.data import DataLoader
from torch.nn import functional as F
from copy import deepcopy
from torch_geometric.nn import GCNConv, global_add_pool, BatchNorm
from tqdm.auto import trange
from active_learning.hyperopt import optimize_hyperparameters


class MLP(torch.nn.Module):
    def __init__(self, in_feats: int = 1024, n_hidden: int = 1024, n_out: int = 2, n_layers: int = 3, seed: int = 42,
                 lr: float = 3e-4, epochs: int = 50, anchored: bool = True, l2_lambda: float = 3e-4,
                 weight_decay: float = 0):
        super().__init__()
        self.seed, self.lr, self.l2_lambda, self.epochs, self.anchored = seed, lr, l2_lambda, epochs, anchored
        self.weight_decay = weight_decay
        torch.manual_seed(seed)

        self.fc = torch.nn.ModuleList()
        self.fc_norms = torch.nn.ModuleList()
        for i in range(n_layers):
            self.fc.append(torch.nn.Linear(in_feats if i == 0 else n_hidden, n_hidden))
            self.fc_norms.append(BatchNorm(n_hidden, allow_single_element=True))
        self.out = torch.nn.Linear(n_hidden, n_out)

    def reset_parameters(self):
        for lin, norm in zip(self.fc, self.fc_norms):
            lin.reset_parameters()
            norm.reset_parameters()
        self.out.reset_parameters()

    def forward(self, x: Tensor) -> Tensor:
        for lin, norm in zip(self.fc, self.fc_norms):
            x = lin(x)
            x = norm(x)
            x = F.relu(x)

        x = self.out(x)
        x = F.log_softmax(x, 1)

        return x


class GCN(torch.nn.Module):
    def __init__(self, in_feats: int = 59, n_hidden: int = 1024, num_conv_layers: int = 3, lr: float = 3e-4,
                 epochs: int = 50, n_out: int = 2, n_layers: int = 3, seed: int = 42, anchored: bool = True,
                 l2_lambda: float = 3e-4, weight_decay: float = 0):

        super().__init__()
        self.seed, self.lr, self.l2_lambda, self.epochs, self.anchored = seed, lr, l2_lambda, epochs, anchored
        self.weight_decay = weight_decay

        self.atom_embedding = torch.nn.Linear(in_feats, n_hidden)

        self.convs = torch.nn.ModuleList()
        self.norms = torch.nn.ModuleList()
        for _ in range(num_conv_layers):
            self.convs.append(GCNConv(n_hidden, n_hidden))
            self.norms.append(BatchNorm(n_hidden, allow_single_element=True))

        self.fc = torch.nn.ModuleList()
        self.fc_norms = torch.nn.ModuleList()
        for i in range(n_layers):
            self.fc.append(torch.nn.Linear(n_hidden, n_hidden))
            self.fc_norms.append(BatchNorm(n_hidden, allow_single_element=True))

        self.out = torch.nn.Linear(n_hidden, n_out)

    def reset_parameters(self):
        self.atom_embedding.reset_parameters()
        for conv, norm in zip(self.convs, self.norms):
            conv.reset_parameters()
            norm.reset_parameters()
        for lin, norm in zip(self.fc, self.fc_norms):
            lin.reset_parameters()
            norm.reset_parameters()
        self.out.reset_parameters()

    def forward(self, x: Tensor, edge_index: Tensor, batch: Tensor) -> Tensor:
        # Atom Embedding:
        x = F.elu(self.atom_embedding(x))

        # Graph convolutions
        for conv, norm in zip(self.convs, self.norms):
            x = conv(x, edge_index)
            x = norm(x)
            x = F.relu(x)

        # Perform global pooling by sum pooling
        x = global_add_pool(x, batch)

        for lin, norm in zip(self.fc, self.fc_norms):
            x = lin(x)
            x = norm(x)
            x = F.relu(x)

        x = self.out(x)
        x = F.log_softmax(x, 1)

        return x


class Model(torch.nn.Module):
    def __init__(self, architecture: str, **kwargs):
        super().__init__()
        assert architecture in ['gcn', 'mlp']
        self.architecture = architecture
        self.model = MLP(**kwargs) if architecture == 'mlp' else GCN(**kwargs)

        self.device_type = "cuda" if torch.cuda.is_available() else "cpu"
        self.device = torch.device(self.device_type)
        self.loss_fn = torch.nn.NLLLoss()

        # Move the whole model to the gpu
        self.model = self.model.to(self.device)
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=self.model.lr,
                                          weight_decay=self.model.weight_decay)

        # Save initial weights in the model for the anchored regularization and move them to the gpu
        if self.model.anchored:
            self.model.anchor_weights = deepcopy({i: j for i, j in self.model.named_parameters()})
            self.model.anchor_weights = {i: j.to(self.device) for i, j in self.model.anchor_weights.items()}

        self.train_loss = []
        self.epochs, self.epoch = self.model.epochs, 0

    def train(self, dataloader: DataLoader, epochs: int = None, verbose: bool = True) -> None:

        bar = trange(self.epochs if epochs is None else epochs, disable=not verbose)
        scaler = torch.cuda.amp.GradScaler()

        for _ in bar:
            running_loss = 0
            items = 0

            for idx, batch in enumerate(dataloader):

                self.optimizer.zero_grad()

                with torch.autocast(device_type=self.device_type, dtype=torch.float16):

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

                    scaler.scale(loss).backward()
                    # loss.backward()
                    scaler.step(self.optimizer)
                    # self.optimizer.step()
                    scaler.update()

                    running_loss += loss.item()
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
            with torch.autocast(device_type=self.device_type, dtype=torch.float16):
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
    def __init__(self, ensemble_size: int = 10, seed: int = 0, architecture: str = 'mlp', **kwargs) -> None:
        self.ensemble_size = ensemble_size
        self.architecture = architecture
        self.seed = seed
        rng = np.random.default_rng(seed=seed)
        self.seeds = rng.integers(0, 1000, ensemble_size)
        self.models = {i: Model(seed=s, architecture=architecture, **kwargs) for i, s in enumerate(self.seeds)}

    def optimize_hyperparameters(self, x, y: DataLoader, **kwargs):
        # raise NotImplementedError
        best_hypers = optimize_hyperparameters(x, y, architecture=self.architecture, **kwargs)
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
