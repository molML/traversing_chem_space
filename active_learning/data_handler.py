
from typing import Union
import numpy as np
import torch
import os
from config import ROOT_DIR


class Handler:
    def __init__(self, n_start: int = 64, bias: str = 'random', seed: int = 42) -> None:

        assert bias in ['random', 'small', 'large'], "'bias' has to be either 'random', 'small', or 'large'"
        assert n_start <= 64, 'Number of starting molecules has to be <= 64'

        self.index_smiles = torch.load(os.path.join(ROOT_DIR, 'data', 'screen', 'index_smiles'))
        self.smiles_index = torch.load(os.path.join(ROOT_DIR, 'data', 'screen', 'smiles_index'))
        self.all_y = torch.load(os.path.join(ROOT_DIR, 'data', 'screen', 'y'))

        self.selected_start_cluster = None
        self.train_idx, self.screen_idx = self.get_start_data(n_start=n_start, bias=bias, seed=seed)
        self.picks = [self.train_idx]

    def get_start_data(self, n_start: int = 64, bias: str = 'random', seed: int = 0) -> (np.ndarray, np.ndarray):

        rng = np.random.default_rng(seed=seed)
        starting_clusters = torch.load(os.path.join(ROOT_DIR, 'data/screen/starting_clusters'))
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





    # def build_train_data(self, n_train: int, balanced: float = False, biased: bool = True, seed: int = 42,
    #                      epsilon: int = 0):
    #
    #     assert not (balanced and biased), "We cannot get a train set that is both balanced and biased"
    #
    #     rng = np.random.default_rng(seed=seed)
    #
    #     if biased:
    #         train_idx = biased_start(n=n_train, epsilon=epsilon, seed=seed)
    #     elif not balanced:
    #         train_idx = rng.integers(0, len(self.all_y), n_train)
    #     else:
    #         # Get all the locations of the actives/inactives
    #         actives_idx = np.where(self.all_y == 1)[0]
    #         inactives_idx = np.where(self.all_y != 1)[0]
    #         # Get random indices of indices to pick
    #         train_actives_idx = rng.integers(0, len(actives_idx), round(n_train * balanced))
    #         train_inactives_idx = rng.integers(0, len(inactives_idx), n_train - round(n_train * balanced))
    #
    #         # Get the actual indices and add them together in one array. Shuffle it
    #         train_idx = np.concatenate((actives_idx[train_actives_idx], inactives_idx[train_inactives_idx]))
    #         rng.shuffle(train_idx)
    #
    #     screen_idx = np.array([i for i in range(len(self.all_y)) if i not in train_idx])
    #     assert len(np.intersect1d(screen_idx, train_idx)) == 0, "Something went wrong selecting train/screen samples"
    #
    #     return train_idx, screen_idx



# def biased_start(n: int = 100, epsilon: float = 0.1, seed: int = 42, name: str = 'screen', root: str = 'data') -> \
#         np.ndarray:
#     """ We start a biased training set with a set of scaffolds based on a hit. This is to simulate a real-case scenario
#     where an optimization series is performed based on a hit compound. We take all the scaffolds from the set of
#     molecules that have same scaffold as the hit. If this set is depleted, we move on towards the next set of scaffolds
#     that is closest (in Tanimoto) to the starting set. To include some different molecules, we also include molecules
#     from random scaffold sets. This is 10% of the total set by default.
#     """
#
#     scaffolds = torch.load(os.path.join(root, name, "scaffolds"))
#     y = torch.load(os.path.join(root, name, "y"))
#     screen_smiles = torch.load(os.path.join(root, name, "smiles"))
#     smiles_index = torch.load(os.path.join(root, name, "smiles_index"))
#     hit_smiles = screen_smiles[np.where(y == 1)]
#
#     rng = np.random.default_rng(seed=seed)
#
#     starting_scaffolds_smiles = []
#     other_scaffolds_smiles = []
#     for scaffold, set_of_smiles in scaffolds.items():
#         if len(set_of_smiles) >= 10 and any([i in hit_smiles for i in set_of_smiles]):
#             starting_scaffolds_smiles.append(list(set_of_smiles))
#         else:
#             other_scaffolds_smiles.append(list(set_of_smiles))
#
#     rng.shuffle(starting_scaffolds_smiles)
#     rng.shuffle(other_scaffolds_smiles)
#
#     # order starting scaffolds on similarity.
#     fps = smiles_to_ecfp([smi_to_scaff(i[0]) for i in starting_scaffolds_smiles], to_array=False)
#     S = BulkTanimotoSimilarity(fps[0], fps)
#     sim_order = np.argsort(S)[::-1]
#     starting_scaffolds_smiles = [starting_scaffolds_smiles[idx] for idx in sim_order]
#
#     # Start the set of selected molecules out with a hit
#     hits_in_set = [i for i in starting_scaffolds_smiles[0] if i in hit_smiles]
#     picked_smiles = {hits_in_set[rng.integers(0, len(hits_in_set))]}
#
#     # add the whole set of scaffolds, if exhausted, go to the next set (they are ranked on similarity to the first set)
#     for scaff_list in starting_scaffolds_smiles:
#         rng.shuffle(scaff_list)
#         for smi in scaff_list:
#             if len(picked_smiles) < n - n*epsilon:
#                 picked_smiles.add(smi)
#
#     # For the remainder, add some other scaffolds (from a random set). Many of them are 1 member sets.
#     for smi_set in other_scaffolds_smiles:
#         rng.shuffle(smi_set)
#         if len(picked_smiles) < n:
#             picked_smiles.add(smi_set[0])
#
#         len(other_scaffolds_smiles[1])
#
#     return np.array([smiles_index[smi] for smi in picked_smiles])
