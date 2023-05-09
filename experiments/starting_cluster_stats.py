

from active_learning.utils import get_tanimoto_matrix
from active_learning.data_prep import MasterDataset
import pandas as pd
from tqdm import tqdm
import numpy as np
from active_learning.data_handler import Handler

if __name__ == '__main__':

    ds_screen = MasterDataset('screen', representation='ecfp')

    starting_samples = {'random_hits': [], 'random_mean_tani': [],
                        'small_hits': [], 'small_mean_tani': [], 'cluster_idx': [],
                        'large_hits': [], 'large_mean_tani': []}

    for s in tqdm(range(1000)):
        for bias in ['random', 'small', 'large']:

            H = Handler(n_start=64, bias=bias, seed=s)
            x, y, smiles = ds_screen[H.get_idx()[0]]
            starting_samples[f"{bias}_hits"].append(sum(y))
            starting_samples[f"{bias}_mean_tani"].append(
                np.mean(get_tanimoto_matrix(smiles, as_vector=True, verbose=False)))
            if bias == 'small':
                starting_samples[f"cluster_idx"].append(H.selected_start_cluster)


    pd.DataFrame(starting_samples).to_csv('starting_sampling.csv')

