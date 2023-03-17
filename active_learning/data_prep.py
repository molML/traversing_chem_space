


# give all molecules an index
# np array with ecfps
# giant torch tensor with all graphs


import pandas as pd
import torch_geometric.data
from tqdm import tqdm
import torch
from torch_geometric.loader import DataLoader
from torch_geometric.data import Data
from utils import molecular_graph_featurizer as smiles_to_graph
from utils import smiles_to_ecfp
import numpy as np
from collections import OrderedDict
from rdkit import Chem
import random


def canonicalize(smiles: str, sanitize: bool = True):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles, sanitize=sanitize))


def get_data(random_state: int = 42):

    # read smiles from file and canonicalize them
    with open('data/inactives.smi') as f:
        inactives = [canonicalize(smi.strip().split()[0]) for smi in f.readlines()]
    with open('data/actives.smi') as f:
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


def split_data(df: pd.DataFrame, random_state: int = 42, screen_test_split: float = 0.9):

    from sklearn.model_selection import train_test_split
    df_screen, df_test = train_test_split(df, stratify=df['y'].tolist(), train_size=screen_test_split,
                                          random_state=random_state)

    # write to csv
    df_screen.to_csv('data/screen_0.csv', index=False)
    df_test.to_csv('data/test_0.csv', index=False)

    return df_screen, df_test


def featurize_data(df_screen, df_test):

    # Compute all molecular graphs and put them in a dict
    screen_graphs = OrderedDict(((s, smiles_to_graph(s, y=y)) for s, y in tqdm(zip(df_screen.smiles, df_screen.y))))
    test_graphs = OrderedDict(((s, smiles_to_graph(s, y=y)) for s, y in tqdm(zip(df_test.smiles, df_test.y))))

    screen_ecfps = OrderedDict([('smiles_index', {smi: i for i, smi in enumerate(df_screen.smiles)}),
                                ('smiles', np.array(df_screen.smiles.tolist())),
                                ('x', smiles_to_ecfp(df_screen.smiles)),
                                ('y', np.array(df_screen.y))])

    test_ecfps = OrderedDict([('smiles_index', {smi: i for i, smi in enumerate(df_test.smiles)}),
                              ('smiles', np.array(df_test.smiles.tolist())),
                              ('x', smiles_to_ecfp(df_test.smiles)),
                              ('y', np.array(df_test.y))])

    torch.save(screen_graphs, 'data/screen_graphs.pt')
    torch.save(test_graphs, 'data/test_graphs.pt')
    torch.save(screen_ecfps, 'data/screen_ecfps.pt')
    torch.save(test_ecfps, 'data/test_ecfps.pt')


if __name__ == '__main__':

    df = get_data()
    df_screen, df_test = split_data(df, screen_test_split=0.9)
    featurize_data(df_screen, df_test)