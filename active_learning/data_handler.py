"""

This script contains the class that manages the data over cycles

    - Handler: class that manages getting the start data and adding new samples to it every active learning cycle.

    Author: Derek van Tilborg, Eindhoven University of Technology, May 2023

"""

from typing import Union
import numpy as np
import torch
import os
from config import ROOT_DIR


class Handler:
    def __init__(self, n_start: int = 64, bias: str = 'random', seed: int = 42, dataset: str = 'ALDH1') -> None:

        assert bias in ['random', 'small', 'large'], "'bias' has to be either 'random', 'small', or 'large'"
        assert n_start <= 64 or bias == 'random', 'Number of starting molecules has to be <= 64'

        self.index_smiles = torch.load(os.path.join(ROOT_DIR, 'data', dataset, 'screen', 'index_smiles'))
        self.smiles_index = torch.load(os.path.join(ROOT_DIR, 'data', dataset, 'screen', 'smiles_index'))
        self.all_y = torch.load(os.path.join(ROOT_DIR, 'data', dataset, 'screen', 'y'))

        self.dataset = dataset
        self.selected_start_cluster = None
        self.train_idx, self.screen_idx = self.get_start_data(n_start=n_start, bias=bias, seed=seed)
        self.picks = [self.train_idx]

    def get_start_data(self, n_start: int = 64, bias: str = 'random', seed: int = 0) -> (np.ndarray, np.ndarray):

        rng = np.random.default_rng(seed=seed)
        starting_clusters = torch.load(os.path.join(ROOT_DIR, f'data/{self.dataset}/screen/starting_clusters'))
        n_clusters = len(starting_clusters)
        self.selected_start_cluster = seed if seed <= n_clusters else rng.integers(0, n_clusters)

        if bias == 'random':
            # get a random hit to start out with
            hits_idx = np.where(self.all_y == 1)[0]
            selected_hit_idx = np.array([hits_idx[rng.integers(0, len(hits_idx))]])

            # get the other random molecules
            remaining_idxs = np.array([i for i in range(len(self.all_y)) if i not in selected_hit_idx])
            selected_others_idx = rng.integers(0, len(remaining_idxs), n_start - 1)
        else:
            # select a random cluster
            cluster_smiles = starting_clusters[self.selected_start_cluster][0 if bias == 'large' else 1]

            # get the molecule indices and labels
            cluster_smiles_idx = np.array([self.smiles_index[smi] for smi in cluster_smiles])
            cluster_smiles_labels = self.all_y[cluster_smiles_idx]

            # get all hits and select a random hit as a starting point
            hits_idx = cluster_smiles_idx[np.where(cluster_smiles_labels == 1)[0]]
            selected_hit_idx = np.array([hits_idx[rng.integers(0, len(hits_idx))]])

            # get the other random molecules from the cluster
            remaining_idxs = np.array([i for i in cluster_smiles_idx if i not in selected_hit_idx])
            selected_others_idx = remaining_idxs[rng.integers(0, len(remaining_idxs), n_start - 1)]

        train_idx = np.concatenate((selected_hit_idx, selected_others_idx))
        rng.shuffle(train_idx)

        screen_idx = np.array([i for i in range(len(self.all_y)) if i not in train_idx])
        assert len(np.intersect1d(screen_idx, train_idx)) == 0, "Something went wrong selecting train/screen samples"

        return train_idx, screen_idx

    def add(self, picks: Union[list, np.ndarray]):
        # Get the corresponding indices of the master dataset, and save it in self.acquired
        added_idx = np.array([self.smiles_index[smi] for smi in picks])
        self.picks.append(added_idx)

        self.train_idx = np.concatenate((self.train_idx, added_idx))
        self.screen_idx = np.array([i for i in range(len(self.all_y)) if i not in self.train_idx])

    def get_idx(self) -> (np.ndarray, np.ndarray):
        return self.train_idx, self.screen_idx

    def __call__(self, *args, **kwargs):
        return self.get_idx()
