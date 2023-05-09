
from torch_geometric.data import Data
from rdkit import Chem
from tqdm import tqdm
from sklearn.metrics import balanced_accuracy_score, roc_auc_score, precision_score, recall_score, confusion_matrix
import numpy as np
import pandas as pd
import torch
from rdkit.Chem import AllChem, DataStructs
from typing import Union, Optional
from torch_geometric.loader import DataLoader as pyg_DataLoader
from torch import Tensor
from torch.utils.data import TensorDataset, DataLoader


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
        self.precision = [0]
        self.tpr = [0]
        self.tn, self.fp, self.fn, self.tp = [0], [0], [0], [0]

    def eval(self, y_hat: torch.Tensor, y: torch.Tensor):

        y = y.cpu() if type(y) is torch.Tensor else torch.tensor(y)
        y_hat = y_hat.cpu() if type(y_hat) is torch.Tensor else torch.tensor(y_hat)

        # calc_binary_accuracy
        acc = torch.sum(y_hat == y) / len(y)
        self.binary_accuracy.append(acc.item())

        # calc_balanced_accuracy
        balanced_acc = balanced_accuracy_score(y, y_hat)
        self.balanced_accuracy.append(balanced_acc)

        # calc_precision
        try:
            self.precision.append(precision_score(y, y_hat, zero_division=0))
        except:
            self.precision.append(0)

        # calc recall
        try:
            self.tpr.append(recall_score(y, y_hat))
        except:
            self.tpr.append(0)

        # calc confusion
        tn, fp, fn, tp = confusion_matrix(y, y_hat).ravel()
        self.tn.append(tn)
        self.fp.append(fp)
        self.fn.append(fn)
        self.tp.append(tp)

    def __repr__(self):
        return f"Binary accuracy:    {self.binary_accuracy[-1]:.4f}\n" \
               f"Balanced accuracy:  {self.balanced_accuracy[-1]:.4f}\n" \
               f"Precision:          {self.precision[-1]:.4f}\n" \
               f"True positive rate: {self.tpr[-1]:.4f}\n" \
               f"Hits:               {self.tp[-1]}\n" \
               f"Misses:             {self.fn[-1]}\n" \
               f"False positives:    {self.fp[-1]}\n" \
               f"True negatives:     {self.tn[-1]}\n"

    def to_dataframe(self, colnames: str = ''):
        df = pd.DataFrame({'cycle': list(range(len(self.tp))), 'binary_accuracy': self.binary_accuracy,
                           'balanced_accuracy': self.balanced_accuracy, 'precision': self.precision, 'tpr': self.tpr,
                           'tp': self.tp, 'fn': self.fn, 'fp': self.fp, 'tn': self.tn})
        df.columns = [f"{colnames}{i}" for i in df.columns]

        return df


def mol_descriptor(smiles: list, scale: bool = True):
    from rdkit.Chem.QED import qed
    from rdkit.Chem import Descriptors, rdMolDescriptors
    from sklearn import preprocessing as pre
    # smiles = smiles[0]
    X = []
    for smi in tqdm(smiles):
        m = Chem.MolFromSmiles(smi)
        x = np.array([Descriptors.TPSA(m),
                      Descriptors.MolLogP(m),
                      Descriptors.MolWt(m),
                      Descriptors.FpDensityMorgan2(m),
                      Descriptors.HeavyAtomMolWt(m),
                      Descriptors.MaxPartialCharge(m),
                      Descriptors.MinPartialCharge(m),
                      Descriptors.NumRadicalElectrons(m),
                      Descriptors.NumValenceElectrons(m),
                      rdMolDescriptors.CalcFractionCSP3(m),
                      rdMolDescriptors.CalcNumRings(m),
                      rdMolDescriptors.CalcNumRotatableBonds(m),
                      rdMolDescriptors.CalcNumLipinskiHBD(m),
                      rdMolDescriptors.CalcNumLipinskiHBA(m),
                      rdMolDescriptors.CalcNumHeterocycles(m),
                      rdMolDescriptors.CalcNumHeavyAtoms(m),
                      rdMolDescriptors.CalcNumAromaticRings(m),
                      rdMolDescriptors.CalcNumAtoms(m),
                      qed(m)])
        X.append(x)

    if scale:
        return pre.MinMaxScaler().fit_transform(np.array(X))
    return np.array(X)


def tsne():

    from active_learning.data_prep import MasterDataset
    from sklearn.manifold import TSNE
    import numpy as np
    import pandas as pd
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.DataStructs import BulkTanimotoSimilarity
    from sklearn import preprocessing as pre

    ds_test = MasterDataset('test', representation='ecfp')
    ds_screen = MasterDataset('screen', representation='ecfp')

    x_test, y_test, smiles_test = ds_test.all()
    x_screen, y_screen, smiles_screen = ds_screen.all()

    x = np.concatenate((x_test, x_screen))
    y = np.concatenate((y_test, y_screen))
    smiles = np.concatenate((smiles_test, smiles_screen))
    split = ['test'] * len(x_test) + ['screen'] * len(x_screen)

    selected = [1 if i in smiles_screen else 0 for i in smiles]

    X_raw = mol_descriptor(smiles, scale=False)
    X_scaled = pre.MinMaxScaler().fit_transform(X_raw)

    x_df = pd.DataFrame(X_raw, columns=["TPSA", "MolLogP", "MolWt", "FpDensityMorgan2", "HeavyAtomMolWt",
                                    "MaxPartialCharge", "MinPartialCharge", "NumRadicalElectrons",
                                    "NumValenceElectrons", "CalcFractionCSP3", "CalcNumRings", "CalcNumRotatableBonds",
                                    "CalcNumLipinskiHBD", "CalcNumLipinskiHBA", "CalcNumHeterocycles",
                                    "CalcNumHeavyAtoms", "CalcNumAromaticRings", "CalcNumAtoms", "qed"])

    # fps = smiles_to_ecfp(smiles, to_array=False)
    # S = np.array([BulkTanimotoSimilarity(fps[i], fps) for i in tqdm(range(len(fps)))])
    # D = 1 - S

    perplexity = 50
    tsne = TSNE(n_components=2, perplexity=perplexity, n_iter=500)
    results = tsne.fit_transform(X_scaled)
    df = pd.DataFrame({"x": results[:, 0], "y": results[:, 1], "label": y, "split": split, 'smiles': smiles, 'selected': selected})
    df = pd.concat([df, x_df], axis=1)
    df.to_csv(f'tsne_perp{perplexity}_500.csv', index=False)


def get_tanimoto_matrix(smiles: list[str], radius: int = 2, nBits: int = 1024, verbose: bool = True,
                        scaffolds: bool = False, zero_diag: bool = True, as_vector: bool = False):
    """ Calculates a matrix of Tanimoto similarity scores for a list of SMILES string"""
    from active_learning.data_prep import smi_to_scaff

    # Make a fingerprint database
    db_fp = {}
    for smi in smiles:
        if scaffolds:
            m = Chem.MolFromSmiles(smi_to_scaff(smi, includeChirality=False))
        else:
            m = Chem.MolFromSmiles(smi)
        fp = AllChem.GetMorganFingerprintAsBitVect(m, radius=radius, nBits=nBits)
        db_fp[smi] = fp

    smi_len = len(smiles)
    m = np.zeros([smi_len, smi_len], dtype=np.float16)  # We use 16-bit floats to prevent giant matrices
    # Calculate upper triangle of matrix
    for i in tqdm(range(smi_len), disable=not verbose):
        for j in range(i, smi_len):
            m[i, j] = DataStructs.TanimotoSimilarity(db_fp[smiles[i]], db_fp[smiles[j]])
    # Fill in the lower triangle without having to loop (saves ~50% of time)
    m = m + m.T - np.diag(np.diag(m))
    # Fill the diagonal with 0's
    if zero_diag:
        np.fill_diagonal(m, 0)
    if as_vector:
        from scipy.spatial.distance import squareform
        m = squareform(m)

    return m


def to_torch_dataloader(x: Union[list, np.ndarray], y: Optional[np.ndarray] = None, **kwargs) -> \
        Union[DataLoader, pyg_DataLoader]:

    if type(x) is np.ndarray:
        assert y is not None, 'No y values provided'
        return DataLoader(TensorDataset(Tensor(x), Tensor(y).unsqueeze(1).type(torch.LongTensor)), **kwargs)
    else:
        return pyg_DataLoader(x, **kwargs)

