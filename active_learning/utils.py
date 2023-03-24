
from torch_geometric.data import Data
from rdkit import Chem
from tqdm import tqdm
from sklearn.metrics import balanced_accuracy_score, roc_auc_score, precision_score, recall_score, confusion_matrix
import numpy as np
import pandas as pd
import torch


def molecular_graph_featurizer(smiles: str, y=None):

    y = torch.tensor([y]).float()

    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    # RDKIT Atom featurization
    xs = []
    for atom in mol.GetAtoms():
        try:
            x = atom_props(atom)
        except:
            # pass
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

    return x


def smiles_to_ecfp(smiles: list[str], radius: int = 2, nbits: int = 1024, silent: bool = True, to_array: bool = True) \
        -> np.ndarray:
    """ Get a Numpy array of ECFPs from a list of SMILES strings """
    from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
    from rdkit.Chem import MolFromSmiles
    from rdkit.DataStructs import ConvertToNumpyArray

    if type(smiles) is str:
        smiles = [smiles]

    fp = [GetMorganFingerprintAsBitVect(MolFromSmiles(s), radius, nBits=nbits) for s in tqdm(smiles, disable=silent)]

    if not to_array:
        return fp

    output = []
    for f in fp:
        arr = np.zeros((1,))
        ConvertToNumpyArray(f, arr)
        output.append(arr)

    return np.asarray(output)


class Evaluate:
    def __init__(self):
        self.binary_accuracy = [0]
        self.balanced_accuracy = [0]
        self.roc_auc = [0]
        self.precision = [0]
        self.tpr = [0]
        self.tn, self.fp, self.fn, self.tp = [0], [0], [0], [0]

    def eval(self, y_hat: torch.Tensor, y: torch.Tensor, threshold: float = 0.5):

        y_hat = torch.tensor(y_hat) if type(y_hat) is np.ndarray else y_hat.cpu()
        y = torch.tensor(y) if type(y) is np.ndarray else y.cpu()

        y_hat_bin = (y_hat > torch.Tensor([threshold])).float() * 1

        # calc_binary_accuracy
        acc = torch.sum(y_hat_bin == y) / len(y)
        self.binary_accuracy.append(acc.item())

        # calc_balanced_accuracy
        balanced_acc = balanced_accuracy_score(y, y_hat_bin)
        self.balanced_accuracy.append(balanced_acc)

        # calc_roc_auc
        try:
            self.roc_auc.append(roc_auc_score(y, y_hat))
        except:
            self.roc_auc.append(0)

        # calc_precision
        try:
            self.precision.append(precision_score(y, y_hat_bin, zero_division=0))
        except:
            self.precision.append(0)

        # calc recall
        try:
            self.tpr.append(recall_score(y, y_hat_bin))
        except:
            self.tpr.append(0)

        # calc confusion
        tn, fp, fn, tp = confusion_matrix(y, y_hat_bin).ravel()
        self.tn.append(tn)
        self.fp.append(fp)
        self.fn.append(fn)
        self.tp.append(tp)

    def __repr__(self):
        return f"Binary accuracy:    {self.binary_accuracy[-1]:.4f}\n" \
               f"Balanced accuracy:  {self.balanced_accuracy[-1]:.4f}\n" \
               f"Precision:          {self.precision[-1]:.4f}\n" \
               f"True positive rate: {self.tpr[-1]:.4f}\n" \
               f"ROC AUC:            {self.roc_auc[-1]:.4f}\n" \
               f"Hits:               {self.tp[-1]}\n" \
               f"Misses:             {self.fn[-1]}\n" \
               f"False positives:    {self.fp[-1]}\n" \
               f"True negatives:     {self.tn[-1]}\n"

    def to_dataframe(self, colnames: str = ''):
        df = pd.DataFrame({'cycle': list(range(len(self.tp))), 'binary_accuracy': self.binary_accuracy,
                           'balanced_accuracy': self.balanced_accuracy, 'precision': self.precision, 'tpr': self.tpr,
                           'roc_auc': self.roc_auc, 'tp': self.tp, 'fn': self.fn, 'fp': self.fp, 'tn': self.tn})
        df.columns = [f"{colnames}{i}" for i in df.columns]

        return df

