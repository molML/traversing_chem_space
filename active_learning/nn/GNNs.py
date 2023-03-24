
import torch
import torch.nn.functional as F
from torch import Tensor
from torch.nn import Linear
from torch_geometric.nn import GCNConv, global_add_pool, GraphNorm
from torch.utils.data import DataLoader
from warnings import warn
from tqdm.auto import trange
from torch_geometric.nn.dense.linear import Linear
from torch_geometric.utils import add_remaining_self_loops, to_dense_adj, degree
from torch_geometric.utils.num_nodes import maybe_num_nodes
from torch_geometric.loader import DataLoader as pyg_DataLoader

# GCNEnsemble
# BayesianGCN
# Every model should have
#
# train(x, y, **kwargs) -> None
# predict(x, **kwargs) -> y_hat, y_hat_sigma, y_hat_mu




class GCN:
    def __init__(self, seed: int = 42, batch_size: int = 128, **kwargs) -> None:

        torch.manual_seed(seed)
        self.model = GCNModel(**kwargs)
        self.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        self.loss_fn = torch.nn.BCELoss()

        # Move the whole model to the gpu
        self.model = self.model.to(self.device)
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=self.model.lr)

        self.train_loss = []
        self.epochs, self.epoch = self.model.epochs, 0
        self.batch_size = batch_size

    def train(self, x, y, epochs: int = None, verbose: bool = True) -> None:

        bar = trange(self.epochs if epochs is None else epochs, disable=not verbose)
        dataloader = pyg_DataLoader(x, batch_size=self.batch_size)
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

    def predict(self, x) -> Tensor:
        """ Predict

        :param dataloader: Torch geometric data loader with data
        :return: A 1D-tensors
        """
        dataloader = pyg_DataLoader(x, batch_size=self.batch_size)
        y_hats = torch.tensor([]).to(self.device)
        with torch.no_grad():
            for batch in dataloader:
                batch.to(self.device)
                y_hat = self.model(batch.x.float(), batch.edge_index, batch.batch)
                y_hats = torch.cat((y_hats, y_hat), 0)

        return y_hats


class GCNModel(torch.nn.Module):
    def __init__(self, in_channels: int = 59, hidden_channels: int = 512, num_conv_layers: int = 3, lr: float = 0.001,
                 epochs: int = 200):
        super(GCNModel, self).__init__()
        self.lr = lr
        self.epochs = epochs
        self.in_channels = in_channels
        self.hidden_channels = hidden_channels
        self.num_conv_layers = num_conv_layers

        self.lin_0 = Linear(in_channels, hidden_channels)
        self.lin_1 = Linear(hidden_channels, hidden_channels)
        self.lin_2 = Linear(hidden_channels, hidden_channels)
        self.lin_out = Linear(hidden_channels, 1)

        self.convs, self.norms = torch.nn.ModuleList(), torch.nn.ModuleList()
        for _ in range(num_conv_layers):
            self.convs.append(GCNConv(hidden_channels, hidden_channels))
            self.norm = GraphNorm(hidden_channels)

    def reset_parameters(self):
        self.lin_0.reset_parameters()
        self.lin_1.reset_parameters()
        self.lin_2.reset_parameters()
        self.lin_out.reset_parameters()
        for conv, norm in zip(self.convs, self.norms):
            conv.reset_parameters()
            norm.reset_parameters()

    def forward(self, x: Tensor, edge_index: Tensor, batch: Tensor) -> Tensor:
        # Atom Embedding:
        x = F.relu(self.lin_0(x))

        # Graph convolutions
        for conv, norm in zip(self.convs, self.norms):
            x = conv(x, edge_index)
            x = norm(x, batch)
            x = F.relu(x)

        # Perform global pooling by sum pooling
        x = global_add_pool(x, batch)

        # Prediction head
        x = F.relu(self.lin_1(x))
        x = F.relu(self.lin_2(x))
        x = self.lin_out(x)

        return torch.sigmoid(x).squeeze()


def gcn_norm(x, edge_index):
    num_nodes = maybe_num_nodes(edge_index)
    A_tilde, _ = add_remaining_self_loops(edge_index, num_nodes=num_nodes)  # A_tilde = A + I
    D_tilde = degree(A_tilde[0]).pow_(-0.5)  # D_tilde = D(A_tilde)
    D_tilde = D_tilde.masked_fill_(D_tilde == float('inf'), 0)  # convert any infs to 0s
    D_tilde = torch.diag(D_tilde)

    return D_tilde @ to_dense_adj(A_tilde).squeeze() @ D_tilde @ x.float()
