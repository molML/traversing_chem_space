
from torch_geometric.data import Data
from rdkit.Chem import rdPartialCharges
from typing import Union
import numpy as np
import sys
import torch
from rdkit import Chem
from tqdm import tqdm


import os.path
from torch_geometric.loader import DataLoader
from sklearn.model_selection import StratifiedKFold
import numpy as np
import pandas as pd
import torch


def bulk_featurizer(smiles: Union[list, str], y: Union[list, float] = None, return_failed: bool = False,
                          silent: bool = False, hide_progress_bar: bool = True):
    if type(smiles) is str:
        smiles = list(str)
        y = list(y)

    graphs, failed_smiles = [], []
    for smi, y_val in tqdm(zip(smiles, y), disable=hide_progress_bar):
        G = molecular_graph_featurizer(smi, y=torch.tensor([y_val]).float())
        if type(G) is str:
            failed_smiles.append(G)
        else:
            graphs.append(G)

    if not silent and len(failed_smiles) > 0:
        print(f"Failed featurizing {len(failed_smiles)}/{len(smiles)} molecules: {failed_smiles}", flush=True,
              file=sys.stderr)

    if return_failed:
        return graphs, failed_smiles

    return graphs


def molecular_graph_featurizer(smiles: str, y=None):

    y = torch.tensor([y]).float()

    mol = Chem.MolFromSmiles(smiles)

    # RDKIT Atom featurization
    xs = []
    for atom in mol.GetAtoms():
        try:
            x = atom_props(atom)
        except:
            return smiles
        xs.append(x)
    x = torch.tensor(xs)

    # Edge featurization
    edge_indices, edge_attrs = [], []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()

        edge_indices += [[i, j], [j, i]]

    edge_index = torch.tensor(edge_indices)
    edge_index = edge_index.t().to(torch.long).view(2, -1)

    # Sort indices.
    if edge_index.numel() > 0:
        perm = (edge_index[0] * x.size(0) + edge_index[1]).argsort()
        edge_index = edge_index[:, perm]

    if torch.isnan(x).any():
        return smiles
        # raise ValueError(f"Featurizing {smiles} gave nan(s)")

    graph = Data(x=x, edge_index=edge_index, smiles=smiles, y=y)

    return graph


def atom_props(atom):

    x = []

    atom_types = ['C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'I', 'P', 'Si', 'B', 'Se']
    symbols = [0] * 12
    symbols[atom_types.index(atom.GetSymbol())] = 1
    x += symbols

    degrees = [0] * 6  # {1, 2, 3, 4, 5, 6}
    degrees[atom.GetDegree()-1] = 1
    x += degrees

    total_degree = [0] * 6  # {1, 2, 3, 4, 5, 6}
    total_degree[atom.GetTotalDegree()-1] = 1
    x += total_degree

    explicit_valance = [0] * 6  # {1, 2, 3, 4, 5, 6}
    explicit_valance[atom.GetExplicitValence()-1] = 1
    x += explicit_valance

    implicit_valence = [0] * 4  # {0, 1, 2, 3}
    implicit_valence[atom.GetImplicitValence()] = 1
    x += implicit_valence

    GetTotalValence = [0] * 6  # {1, 2, 3, 4, 5, 6}
    GetTotalValence[atom.GetImplicitValence()-1] = 1
    x += GetTotalValence

    implicit_Hs = [0] * 4  # {0, 1, 2, 3}
    implicit_Hs[atom.GetNumImplicitHs()] = 1
    x += implicit_Hs

    total_Hs = [0] * 4  # {0, 1, 2, 3}
    total_Hs[atom.GetTotalNumHs()] = 1
    x += total_Hs

    formal_charge = [0] * 5  # {-1, 0, 1, 2, 3}
    formal_charge[atom.GetFormalCharge()+1] = 1
    x += formal_charge

    hybridization = [0] * 6
    possible_hybridizations = ['SP', 'SP2', 'SP3', 'SP2D', 'SP3D', 'SP3D2']
    hybridization[possible_hybridizations.index(atom.GetHybridization().name)] = 1
    x += hybridization

    partial_charge = [sigmoid(float(atom.GetProp('_GasteigerCharge')))]
    x.append(partial_charge)

    atomic_weight = sigmoid(Chem.GetPeriodicTable().GetAtomicWeight(atom.GetSymbol()))
    x.append(atomic_weight)

    return x


def sigmoid(x: float):
    return 1.0 / (1.0 + float(np.exp(-x)))


def molecule_from_smiles(smiles: str):
    """ Sanitize a molecule from a SMILES string"""

    molecule = Chem.MolFromSmiles(smiles, sanitize=False)

    # If sanitization is unsuccessful, catch the error, and try again without
    # the sanitization step that caused the error
    flag = Chem.SanitizeMol(molecule, catchErrors=True)
    if flag != Chem.SanitizeFlags.SANITIZE_NONE:
        Chem.SanitizeMol(molecule, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ flag)

    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    Chem.rdPartialCharges.ComputeGasteigerCharges(molecule)

    return molecule

# dataset_path = 'data/screen.csv'
# def load_data(dataset_path, batch_size: int = 32, shuffle: bool = False):
#
#     pickled = dataset_path.replace('.csv', '.pt')
#     if os.path.exists(pickled):
#         data = torch.load(pickled)
#     else:
#         df = pd.read_csv(dataset_path)
#         graphs = bulk_featurizer(df.smiles.tolist(), df.y.tolist(), hide_progress_bar=False)
#
#         data = DataLoader(graphs, batch_size=batch_size, shuffle=shuffle)
#         torch.save(data, pickled)
#
#     return data

# indices_path = 'data/split_indices.pt'
def load_data(indices_path: str = 'data/split_indices.pt', data_path: str = 'data/molecular_graphs.pt',
              batch_size: int = 32, type: str = 'graph'):
    split_indices = torch.load(indices_path)

    if type == 'graph':
        data = torch.load(data_path)

        screen = DataLoader([data[idx] for idx in split_indices['screen']], batch_size=batch_size, shuffle=False)
        train = DataLoader([data[idx] for idx in split_indices['train']], batch_size=batch_size, shuffle=True)
        test = DataLoader([data[idx] for idx in split_indices['test']], batch_size=batch_size, shuffle=False)
        val = DataLoader([data[idx] for idx in split_indices['val']], batch_size=batch_size, shuffle=False)

        return screen, train, test, val

    if type == 'ECFP':
        data = pd.read_csv('data/all.csv')

        [data[idx] for idx in split_indices['test']]


def smiles_to_ecfp(smiles: list[str], radius: int = 2, nbits: int = 1024, silent: bool = True) -> np.ndarray:
    """ Get a Numpy array of ECFPs from a list of SMILES strings """
    from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
    from rdkit.Chem import MolFromSmiles
    from rdkit.DataStructs import ConvertToNumpyArray

    if type(smiles) is str:
        smiles = [smiles]

    fp = [GetMorganFingerprintAsBitVect(MolFromSmiles(s), radius, nBits=nbits) for s in tqdm(smiles, disable=silent)]

    output = []
    for f in fp:
        arr = np.zeros((1,))
        ConvertToNumpyArray(f, arr)
        output.append(arr)

    return np.asarray(output)


class Evaluate:
    def __init__(self, y_hat: torch.Tensor, y: torch.Tensor):
        self.y_hat = y_hat.to('cpu')
        self.y = y.to('cpu')

        self.binary_accuracy = None
        self.balanced_accuracy = None
        self.roc_auc = None
        self.precision = None
        self.tpr = None
        self.tn, self.fp, self.fn, self.tp = None, None, None, None

    def calc_binary_accuracy(self, threshold: float = 0.5):
        y_hat = (self.y_hat > torch.Tensor([threshold])).float() * 1

        acc = torch.sum(y_hat == self.y) / len(self.y)

        self.binary_accuracy = acc.item()
        return self.binary_accuracy

    def calc_balanced_accuracy(self, threshold: float = 0.5):
        from sklearn.metrics import balanced_accuracy_score

        y_hat = (self.y_hat > torch.Tensor([threshold])).float() * 1
        balanced_acc = balanced_accuracy_score(self.y, y_hat)

        self.balanced_accuracy = balanced_acc
        return self.balanced_accuracy

    def calc_roc_auc(self):
        from sklearn.metrics import roc_auc_score

        self.roc_auc = roc_auc_score(self.y, self.y_hat)
        return self.roc_auc

    def calc_precision(self, threshold: float = 0.5):
        from sklearn.metrics import precision_score

        y_hat = (self.y_hat > torch.Tensor([threshold])).float() * 1
        self.precision = precision_score(self.y, y_hat)

        return self.precision

    def calc_recall(self, threshold: float = 0.5):
        from sklearn.metrics import recall_score

        y_hat = (self.y_hat > torch.Tensor([threshold])).float() * 1
        self.tpr = recall_score(self.y, y_hat)

        return self.tpr

    def calc_confusion(self, threshold: float = 0.5):
        from sklearn.metrics import confusion_matrix

        y_hat = (self.y_hat > torch.Tensor([threshold])).float() * 1
        self.tn, self.fp, self.fn, self.tp = confusion_matrix(self.y, y_hat).ravel()

    def calc_hits(self, **kwargs):
        if self.tp is None:
            self.calc_confusion(**kwargs)
        return self.tp

    def calc_misses(self, **kwargs):
        if self.fn is None:
            self.calc_confusion(**kwargs)
        return self.fn

    def eval(self, silent: bool = False, threshold: float = 0.5):
        self.calc_binary_accuracy(threshold=threshold)
        self.calc_balanced_accuracy(threshold=threshold)
        self.calc_precision(threshold=threshold)
        self.calc_recall(threshold=threshold)
        self.calc_roc_auc()
        self.calc_hits(threshold=threshold)
        self.calc_misses(threshold=threshold)

        if not silent:
            print(self)

    def __repr__(self):
        return f"Binary accuracy:    {self.binary_accuracy:.4f}\n" \
               f"Balanced accuracy:  {self.balanced_accuracy:.4f}\n" \
               f"Precision:          {self.precision:.4f}\n" \
               f"True positive rate: {self.tpr:.4f}\n" \
               f"ROC AUC:            {self.roc_auc:.4f}\n" \
               f"Hits:               {self.tp}\n" \
               f"Misses:             {self.fn}\n" \
               f"False positives:    {self.fp}\n" \
               f"True negatives:     {self.tn}"

    def metrics_fstring(self):
        if self.binary_accuracy is None:
            self.eval(silent=True)
        return f"{self.binary_accuracy},{self.balanced_accuracy},{self.precision},{self.tpr},{self.roc_auc}," \
               f"{self.tp},{self.fn},{self.fp},{self.tn}"
