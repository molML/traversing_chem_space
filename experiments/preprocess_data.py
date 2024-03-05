
import os
import warnings
from collections import Counter
import torch
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from active_learning.data_prep import MasterDataset, load_hdf5, get_data, split_data, similarity_vectors
from config import ROOT_DIR
warnings.simplefilter(action='ignore', category=FutureWarning)


if __name__ == '__main__':

    # Process the data
    for dataset in ['ALDH1', 'PKM2', 'VDR']:

        df = get_data(dataset=dataset)
        df_screen, df_test = split_data(df, screen_size=1000, test_size=200, dataset=dataset)

        MasterDataset(name='screen', df=df_screen, overwrite=True, dataset=dataset)
        MasterDataset(name='test', df=df_test, overwrite=True, dataset=dataset)

        df_screen = pd.read_csv(os.path.join(ROOT_DIR, f'data/{dataset}/original/screen.csv'))
        df_test = pd.read_csv(os.path.join(ROOT_DIR, f'data/{dataset}/original/test.csv'))

        similarity_vectors(df_screen, df_test, dataset=dataset)

    # Perform clustering for each dataset
    for dataset, tani_cutoffs in zip(['PKM2', 'VDR', 'ALDH1'], [[0.80, 0.61], [0.80, 0.70], [0.80, 0.60]]):
        ds_screen = MasterDataset('screen', representation='ecfp', dataset=dataset)
        x_screen, y_screen, smiles_screen = ds_screen.all()
        smiles_index = torch.load(f'data/{dataset}/screen/smiles_index')
        min_supercluster_size = 128
        min_subcluster_size = 64
        n = len(smiles_screen)
        subcluster_mu = np.mean(y_screen.tolist()) * min_subcluster_size
        subcluster_sigma = np.std(y_screen.tolist()) * np.sqrt(min_subcluster_size)

        D = load_hdf5(f'data/{dataset}/screen/tanimoto_distance_vector')

        # Average clustering -> works probably slightly better than complete, as it gives much larger clusters
        linkage = hierarchy.average(D)
        del D

        # Cut the tree in clusters where the average intra-cluster Tanimoto distance is 0.8, 0.6
        cut_clusters = hierarchy.cut_tree(linkage, height=tani_cutoffs)
        cut = np.concatenate((np.array([range(n)]).T, cut_clusters), axis=1)

        # Find the big superclusters
        super_clusters = [clust for clust, cnt in Counter(cut[:, 1]).items() if cnt >= min_supercluster_size]
        cut = cut[[True if i in super_clusters else False for i in cut[:, 1]]]

        # find the subclusters
        sub_clusters = [clust for clust, cnt in Counter(cut[:, -1]).items() if cnt >= min_subcluster_size]


        # put the subclusters and superclusters together
        cluster_smiles = []
        for sub_clust in sub_clusters:
            super_clust = cut[:, 1][cut[:, 2] == sub_clust][0]
            cluster_smiles.append([
                smiles_screen[cut[:, 0][np.where(cut[:, 2] == sub_clust)]],
                smiles_screen[cut[:, 0][np.where(cut[:, 1] == super_clust)]]
            ])
        cluster_smiles = np.array(cluster_smiles, dtype=object)  # len = 46  len(cluster_smiles)

        # Only keep the clusters where the subcluster actually contains a hit
        cluster_smiles_with_hits = []
        for subset, superset in cluster_smiles:
            y_subset = y_screen[np.array([smiles_index[smi] for smi in subset])]
            # print(sum(y_subset))
            if sum(y_subset) > subcluster_mu + subcluster_sigma:
                cluster_smiles_with_hits.append([subset, superset])
        cluster_smiles_with_hits = np.array(cluster_smiles_with_hits, dtype=object)   # len(cluster_smiles_with_hits)

        for i in cluster_smiles_with_hits:
            print(len(i[0]), len(i[1]), len(i[0]) / len(i[1]))

        only_child = []
        for i in range(len(cluster_smiles_with_hits)):
            supercluster = cluster_smiles_with_hits[i][1]
            subcluster = cluster_smiles_with_hits[i][0]

            contains = 0
            if len(np.intersect1d(supercluster, subcluster)) > 0:
                contains = 1

            for j in range(len(cluster_smiles_with_hits)):
                if i != j and len(np.intersect1d(cluster_smiles_with_hits[j][1], subcluster)) > 0:
                    contains += 1

            only_child.append(contains)

        cluster_smiles_with_hits = cluster_smiles_with_hits[np.where(np.array(only_child) == 1)]
        print(len(cluster_smiles_with_hits))

        torch.save(cluster_smiles_with_hits, f'data/{dataset}/screen/starting_clusters')
