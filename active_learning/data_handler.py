

import numpy as np
import torch
import os
from typing import Union


class Handler:
    def __init__(self, n_train: int = 1000, balanced: float = False, seed: int = 42):

        self.index_smiles = torch.load(os.path.join('data', 'screen', 'index_smiles'))
        self.smiles_index = torch.load(os.path.join('data', 'screen', 'smiles_index'))
        self.all_y = torch.load(os.path.join('data', 'screen', 'y'))

        self.train_idx, self.screen_idx = self.build_train_data(n_train=n_train, balanced=balanced, seed=seed)
        self.picks = [self.train_idx]

    def build_train_data(self, n_train: int, balanced: float = False, seed: int = 42):

        rng = np.random.default_rng(seed=seed)

        if not balanced:
            train_idx = rng.integers(0, len(self.all_y), n_train)
        else:
            # Get all the locations of the actives/inactives
            actives_idx = np.where(self.all_y == 1)[0]
            inactives_idx = np.where(self.all_y != 1)[0]
            # Get random indices of indices to pick
            train_actives_idx = rng.integers(0, len(actives_idx), round(n_train * balanced))
            train_inactives_idx = rng.integers(0, len(inactives_idx), n_train - round(n_train * balanced))

            # Get the actual indices and add them together in one array. Shuffle it
            train_idx = np.concatenate((actives_idx[train_actives_idx], inactives_idx[train_inactives_idx]))
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
