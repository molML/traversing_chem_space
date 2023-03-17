

from active_learning.utils import molecular_graph_featurizer as smiles_to_graph
from active_learning.utils import smiles_to_ecfp
import pandas as pd
import numpy as np
import torch
import os
import sys
from collections import OrderedDict
from tqdm import tqdm
from rdkit import Chem


def canonicalize(smiles: str, sanitize: bool = True):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles, sanitize=sanitize))


def get_data(random_state: int = 42):

    # read smiles from file and canonicalize them
    with open('data/original/inactives.smi') as f:
        inactives = [canonicalize(smi.strip().split()[0]) for smi in f.readlines()]
    with open('data/original/actives.smi') as f:
        actives = [canonicalize(smi.strip().split()[0]) for smi in f.readlines()]

    # remove duplicates:
    inactives = list(set(inactives))
    actives = list(set(actives))

    # remove intersecting molecules:
    intersecting_mols = np.intersect1d(inactives, actives)
    inactives = [smi for smi in inactives if smi not in intersecting_mols]
    actives = [smi for smi in actives if smi not in intersecting_mols]

    # add to df
    df = pd.DataFrame({'smiles': inactives + actives,
                       'y': [0] * len(inactives) + [1] * len(actives)})

    # shuffle
    df = df.sample(frac=1, random_state=random_state).reset_index(drop=True)

    return df


def split_data(df: pd.DataFrame, random_state: int = 42, screen_size: float = 0.9) -> (pd.DataFrame, pd.DataFrame):

    from sklearn.model_selection import train_test_split
    df_screen, df_test = train_test_split(df, stratify=df['y'].tolist(), train_size=screen_size,
                                          random_state=random_state)

    # write to csv
    df_screen.to_csv('data/original/screen.csv', index=False)
    df_test.to_csv('data/original/test.csv', index=False)

    return df_screen, df_test


class dataset:
    """ Dataset that holds all data in an indexable way """
    def __init__(self, df: pd.DataFrame, name: str, representation: str = 'ecfp', root: str = 'data') -> None:

        assert representation in ['ecfp', 'graph'], f"'representation' must be 'ecfp' or 'graph', not {representation}"
        self.representation = representation
        self.pth = os.path.join(root, name)

        # If not done already, process all data. Else just load it
        if not os.path.exists(self.pth):
            os.makedirs(os.path.join(root, name))
            self.process(df)
            self.smiles_index, self.index_smiles, self.smiles, self.x, self.y, self.graphs = self.load()
        else:
            self.smiles_index, self.index_smiles, self.smiles, self.x, self.y, self.graphs = self.load()

    def process(self, df: pd.DataFrame) -> None:

        print('Processing data ... ', flush=True, file=sys.stderr)

        index_smiles = OrderedDict({i: smi for i, smi in enumerate(df.smiles)})
        smiles_index = OrderedDict({smi: i for i, smi in enumerate(df.smiles)})
        smiles = np.array(df.smiles)
        x = smiles_to_ecfp(df.smiles, silent=False)
        y = np.array(df.y)
        graphs = [smiles_to_graph(smi, y=y) for smi, y in tqdm(zip(df.smiles, df.y))]

        torch.save(index_smiles, os.path.join(self.pth, 'index_smiles'))
        torch.save(smiles_index, os.path.join(self.pth, 'smiles_index'))
        torch.save(smiles, os.path.join(self.pth, 'smiles'))
        torch.save(x, os.path.join(self.pth, 'x'))
        torch.save(y, os.path.join(self.pth, 'y'))
        torch.save(graphs, os.path.join(self.pth, 'graphs'))

    def load(self) -> (dict, dict, np.ndarray, np.ndarray, np.ndarray, list):

        print('Loading data ... ', flush=True, file=sys.stderr)

        index_smiles = torch.load(os.path.join(self.pth, 'index_smiles'))
        smiles_index = torch.load(os.path.join(self.pth, 'smiles_index'))
        smiles = torch.load(os.path.join(self.pth, 'smiles'))
        x = torch.load(os.path.join(self.pth, 'x'))
        y = torch.load(os.path.join(self.pth, 'y'))
        graphs = torch.load(os.path.join(self.pth, 'graphs'))

        return smiles_index, index_smiles, smiles, x, y, graphs

    def __len__(self) -> int:
        return len(self.smiles)

    def __getitem__(self, idx):
        if type(idx) is int:
            idx = [idx]
        if self.representation == 'ecfp':
            return self.x[idx], self.y[idx]
        if self.representation == 'graph':
            return [self.graphs[i] for i in idx], self.y[idx]


if __name__ == '__main__':

    df = get_data()
    df_screen, df_test = split_data(df, screen_size=0.9)

    dataset(df_screen, 'screen')
    dataset(df_test, 'test', 'graph')
