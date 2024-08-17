import pandas as pd
import numpy as np
from active_learning.nn import Ensemble
from active_learning.data_prep import MasterDataset
from active_learning.data_handler import Handler
from active_learning.utils import Evaluate, to_torch_dataloader
from active_learning.acquisition import Acquisition, logits_to_pred
from tqdm.auto import tqdm
import torch
from torch.utils.data import WeightedRandomSampler
from math import ceil
from rdkit.DataStructs import BulkTanimotoSimilarity
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.Chem import Descriptors
from warnings import warn
from rdkit.DataStructs import ConvertToNumpyArray
from rdkit import Chem
import seaborn as sns
import matplotlib.pyplot as plt


def rdkit_to_array(fp: list) -> np.ndarray:
    """ Convert a list of RDkit fingerprint objects into a numpy array """
    output = []
    for f in fp:
        arr = np.zeros((1,))
        ConvertToNumpyArray(f, arr)
        output.append(arr)
    return np.asarray(output)


def mols_to_ecfp(mols: list, radius: int = 2, nbits: int = 1024, progressbar: bool = False,
                 to_array: bool = False):
    """ Get ECFPs from a list of RDKit molecule objects

    :param mols: list of RDKit mol objects, e.g., as obtained through smiles_to_mols()
    :param radius: Radius of the ECFP (default = 2)
    :param nbits: Number of bits (default = 1024)
    :param progressbar: toggles progressbar (default = False)
    :param to_array: Toggles conversion of RDKit fingerprint objects to a Numpy Array (default = False)
    :return: list of RDKit ECFP fingerprint objects, or a Numpy Array of ECFPs if to_array=True
    """
    fp = [GetMorganFingerprintAsBitVect(m, radius, nBits=nbits) for m in tqdm(mols, disable=not progressbar)]
    if not to_array:
        return fp
    return rdkit_to_array(fp)



mol_sim = []
kind = []
dataset = []
for dataset_name in ['ALDH1', 'PKM2', 'VDR']:
    ds_screen = MasterDataset('screen', representation='ecfp', dataset=dataset_name)

    actives = ds_screen.smiles[ds_screen.y == 1]
    inactives = ds_screen.smiles[ds_screen.y == 0]
    actives_fps = mols_to_ecfp([Chem.MolFromSmiles(smi) for smi in actives], to_array=False)
    inactives_fps = mols_to_ecfp([Chem.MolFromSmiles(smi) for smi in inactives], to_array=False)


    for fp_i in tqdm(range(len(actives_fps))):
        mol_sim.append(np.mean(BulkTanimotoSimilarity(actives_fps[fp_i], actives_fps[fp_i + 1:])))
        kind.append('active - active')
        dataset.append(dataset_name)

    for fp_i in tqdm(range(len(actives_fps))):
        mol_sim.append(np.mean(BulkTanimotoSimilarity(actives_fps[fp_i], inactives_fps)))
        kind.append('active - inactive')
        dataset.append(dataset_name)

df = pd.DataFrame({'Tanimoto similarity': mol_sim, 'kind': kind, 'dataset': dataset})
df.to_csv('figures/data/similarity_histograms.csv', index=False)


data_wide = df.pivot(columns = 'kind',
                     values = 'Tanimoto similarity')

data_wide.plot.kde(figsize = (8, 6),  linewidth = 4)

plt.show()
