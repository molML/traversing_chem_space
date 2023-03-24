

import torch
from torch import nn
import torch.nn.functional as F
from torch import Tensor, tensor
from torch_geometric.nn import GCNConv, global_add_pool, GraphNorm
import pyro
import pyro.distributions as dist
from pyro.nn import PyroModule, PyroSample
from active_learning.nn.GNNs import gcn_norm


class BNN(PyroModule):
    """ Simple feedforward Bayesian NN. We use PyroSample to place a prior over the weight and bias parameters,
     instead of treating them as fixed learnable parameters. See http://pyro.ai/examples/bayesian_regression.html"""

    def __init__(self, in_feats=1024, hidden_size: int = 256, lr: float = 0.001, epochs: int = 10000,
                 to_gpu: bool = True) -> None:
        super().__init__()
        # Define some vars
        self.lr = lr
        self.epochs = epochs
        self.to_gpu = to_gpu
        nh0, nh, nhout = in_feats, hidden_size, 2
        self.device = torch.device("cuda:0" if torch.cuda.is_available() and to_gpu else "cpu")

        self.fc0 = PyroModule[nn.Linear](nh0, nh)
        self.fc0.weight = PyroSample(dist.Normal(0., tensor(1.0, device=self.device)).expand([nh, nh0]).to_event(2))
        self.fc0.bias = PyroSample(dist.Normal(0., tensor(1.0, device=self.device)).expand([nh]).to_event(1))

        self.fc1 = PyroModule[nn.Linear](nh, nh)
        self.fc1.weight = PyroSample(dist.Normal(0., tensor(1.0, device=self.device)).expand([nh, nh]).to_event(2))
        self.fc1.bias = PyroSample(dist.Normal(0., tensor(1.0, device=self.device)).expand([nh]).to_event(1))

        self.fc2 = PyroModule[nn.Linear](nh, nh)
        self.fc2.weight = PyroSample(dist.Normal(0., tensor(1.0, device=self.device)).expand([nh, nh]).to_event(2))
        self.fc2.bias = PyroSample(dist.Normal(0., tensor(1.0, device=self.device)).expand([nh]).to_event(1))

        self.out = PyroModule[nn.Linear](nh, nhout)
        self.out.weight = PyroSample(dist.Normal(0., tensor(1.0, device=self.device)).expand([nhout, nh]).to_event(2))
        self.out.bias = PyroSample(dist.Normal(0., tensor(1.0, device=self.device)).expand([nhout]).to_event(1))

    def forward(self, x: tensor, y: tensor = None) -> tensor:
        # Predict the logits
        x = F.relu(self.fc0(x))
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        logits = self.out(x)

        with pyro.plate("data", x.shape[0]):
            obs = pyro.sample("obs", dist.Categorical(logits=logits), obs=y)

        return logits


class BayesianGCNModel(PyroModule):
    def __init__(self, in_feats=59, hidden_size: int = 128, to_gpu: bool = True) -> None:
        super().__init__()
        # Define some vars
        nh0, nh, nhout = in_feats, hidden_size, 1
        self.device = torch.device("cuda:0" if torch.cuda.is_available() and to_gpu else "cpu")

        self.fc0 = PyroModule[nn.Linear](nh0, nh)
        self.fc0.weight = PyroSample(dist.Normal(0., tensor(1.0, device=self.device)).expand([nh, nh0]).to_event(2))
        self.fc0.bias = PyroSample(dist.Normal(0., tensor(1.0, device=self.device)).expand([nh]).to_event(1))

        self.gcn0 = PyroModule[nn.Linear](nh, nh)
        self.gcn0.weight = PyroSample(dist.Normal(0., tensor(1.0, device=self.device)).expand([nh, nh]).to_event(2))
        self.gcn0.bias = PyroSample(dist.Normal(0., tensor(1.0, device=self.device)).expand([nh]).to_event(1))

        self.gcn1 = PyroModule[nn.Linear](nh, nh)
        self.gcn1.weight = PyroSample(dist.Normal(0., tensor(1.0, device=self.device)).expand([nh, nh]).to_event(2))
        self.gcn1.bias = PyroSample(dist.Normal(0., tensor(1.0, device=self.device)).expand([nh]).to_event(1))

        self.fc1 = PyroModule[nn.Linear](nh, nh)
        self.fc1.weight = PyroSample(dist.Normal(0., tensor(1.0, device=self.device)).expand([nh, nh]).to_event(2))
        self.fc1.bias = PyroSample(dist.Normal(0., tensor(1.0, device=self.device)).expand([nh]).to_event(1))

        self.out = PyroModule[nn.Linear](nh, nhout)
        self.out.weight = PyroSample(dist.Normal(0., tensor(1.0, device=self.device)).expand([nhout, nh]).to_event(2))
        self.out.bias = PyroSample(dist.Normal(0., tensor(1.0, device=self.device)).expand([nhout]).to_event(1))

    def forward(self, x: Tensor, edge_index: Tensor, batch: Tensor, y: Tensor = None) -> tensor:
        # Predict the logits
        x = F.relu(self.fc0(x))
        x = F.relu(self.gcn0(gcn_norm(x, edge_index)))
        x = F.relu(self.gcn1(gcn_norm(x, edge_index)))
        x = global_add_pool(x, batch)
        x = F.relu(self.fc1(x))
        logits = self.out(x).squeeze()

        with pyro.plate("data", logits.shape[0]):
            obs = pyro.sample("obs", dist.Categorical(logits=logits), obs=y)

        return logits


#
# from active_learning.data_prep import MasterDataset
# #
# ds_test = MasterDataset('test', representation='graph')
# x_train, y_train, smiles_train = ds_test[range(100)]
#
# x = pyg_DataLoader(x_train, batch_size=32)
#
# # init stuff
# model = BayesianGCNModel(to_gpu=False)
# guide = AutoDiagonalNormal(model)
# adam = pyro.optim.Adam({"lr": 0.001})
# svi = SVI(model, guide, adam, loss=Trace_ELBO())
#
# # train
# pyro.clear_param_store()
#
# for epoch in range(1000):
#     running_loss = 0.0
#     n_samples = 0
#     for batch in x:
#         # break
#         # svi.step(batch)
#         # x = batch[0]
#         # y = batch[1]
#
#         # ELBO gradient and add loss to running loss
#         running_loss += svi.step(batch)
#         n_samples += len(batch)
#
#     loss = running_loss / n_samples
#     print(loss)


