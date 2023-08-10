
import torch
from active_learning.data_prep import MasterDataset, load_hdf5
import numpy as np
from scipy.cluster import hierarchy
from collections import Counter
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

if __name__ == '__main__':

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
        # linkage = hierarchy.average(D)
        # del D

        linkage = torch.load(f'data/{dataset}/screen/average_linkage_clustering')

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

    # only child clusters [0, 2, 4, 6, 7, 8, 9, 10, 13, 15]
    # new ordering        [0, 1, 2, 3, 4, 5, 6, 7,  8,  9]

acquisitions = ['random', 'exploration', 'exploitation', 'bald', 'similarity', 'exploitation_nolearn']  #  'dynamic'
biases = ["random", "small", "large"]
batch_sizes = [64, 32, 16]
architectures = ["gcn", "mlp"]
retrains = ['True', 'False']

jobs = []
for dataset in ['PKM2', 'VDR']:
    for acq_size in [64]:
        for acq_func in ['random', 'exploration', 'exploitation', 'bald', 'similarity', 'exploitation_nolearn']:
            for arch in ["gcn", "mlp"]:
                for bias in ["random", "small", "large"]:
                    jobs.append(f"{dataset}_{acq_size}_{acq_func}_{arch}_{bias}")

len(jobs)

# # Ward clustering of scaffold similarity
#
# n_super_clusters = 10
# n_sub_clusters = 10
#
# linkage = hierarchy.ward(D)
# clusters = hierarchy.cut_tree(linkage, n_clusters=[n_super_clusters])
# clusters = np.concatenate((np.array([range(n)]).T, clusters), axis=1)
#
# all_sub_clusters = np.zeros((n, 1), dtype=int)
# for i in range(n_super_clusters):
#     smiles_idx = clusters[:, 0][np.where(clusters[:, 1] == i)]
#     smiles = smiles_screen[smiles_idx]
#     D = get_tanimoto_matrix(smiles, as_vector=True, scaffolds=True)
#     Z = hierarchy.ward(D)
#     sub_clusters = hierarchy.cut_tree(Z, n_clusters=[n_sub_clusters])
#     all_sub_clusters[smiles_idx] = sub_clusters
#
#     print(Counter(sub_clusters[:, 0]))
#
# clusters = np.concatenate((clusters, all_sub_clusters), axis=1)
#
#
# smiles_per_cluster = [smiles_screen[cut[:, 0][np.where(cut[:, 1] == i)]] for i in range(10)]
#
#
# for smiles in smiles_per_cluster:
#     D = get_tanimoto_matrix(smiles, as_vector=True, scaffolds=True)
#     Z = hierarchy.ward(D)
#     sub_clusters = hierarchy.cut_tree(Z, n_clusters=[10])
#
#
# # Complete agglomarative clutering with two distance cutoffs. The strict cutoff is a subset of the relaxed cutoff.
#
# linkage = hierarchy.complete(D)
# del D
#
# cut_clusters = hierarchy.cut_tree(linkage, height=[0.8, 0.2])    # [sim_leaves]
# cut = np.concatenate((np.array([range(n)]).T, cut_clusters), axis=1)
#
# big_cluster_size = 100
# big_clusters = [clust for clust, cnt in Counter(cut[:, 1]).items() if cnt >= big_cluster_size]
#
# cut = cut[[True if i in big_clusters else False for i in cut[:, 1]]]
#
# big_stric_clusters = [clust for clust, cnt in Counter(cut[:, -1]).items() if cnt >= 50]
# cut_strict = cut[[True if i in big_stric_clusters else False for i in cut[:, -1]]]
#
# cluster_smiles = []
# for sub_clust in big_stric_clusters:
#     super_clust = cut[:, 1][cut[:, 2] == sub_clust][0]
#
#     cluster_smiles.append([
#         smiles_screen[cut[:, 0][np.where(cut[:, 2] == sub_clust)]],
#         smiles_screen[cut[:, 0][np.where(cut[:, 1] == super_clust)]]
#     ])
#
# for i in cluster_smiles:
#     print(len(i[0]), len(i[1]), len(i[0]) / len(i[1]))
#
# print(sum([1 for i in cluster_smiles if len(i[0]) / len(i[1]) <= 0.5]), len(cluster_smiles))
#
#
#
#
#
#
# cluster_smiles[-1][0]
# cluster_smiles[-1][-1]
#
#
# cut[:, 0][cut[:, 2] == 1]
#
# 28 ;161
#
# pd.DataFrame(cut, columns=['idx', 'clust_08', 'clust_06']).to_csv('tree_clusts.csv', index=False)
#
# len(big_stric_clusters)
#
# Counter(cut[:, 3])
#
# Counter(cut[:, 1])
#
# max(cut[:, 1])
#
# np.array([range(n)]).T.shape
#
#
# smiles_screen[6]
#
#
# cut_clusters = hierarchy.cut_tree(linkage, height=[median_dist])
#
# # hierarchy.dendrogram(linkage, count_sort=True, color_threshold=median_dist, labels=cut_clusters[:, -1], leaf_rotation=0)
# # plt.show()
# distance_cutoff = 0.4
#
#
# cut = np.array((np.array(range(n)), cut_clusters[:, -1], )).T
# nodes_dist_threshold = linkage[linkage[:, 2] < distance_cutoff][:, :2].flatten()
# leaves_dist_threshold = nodes_dist_threshold[nodes_dist_threshold < n]
#
# leaves_cluster = cut[[True if i in leaves_dist_threshold else False for i in cut[:, 0]]]
# # rebase the cluster labels
# for i, c in enumerate(set(leaves_cluster[:, 1])):
#     leaves_cluster[:, 1][leaves_cluster[:, 1] == c] = i
#
# big_clusters = [k for k, v in Counter(leaves_cluster[:, 1]).items() if v > big_cluster_size]
#
# # lif of smiles for each big clusters
# smiles_per_cluster = [smiles_screen[leaves_cluster[:, 0][np.where(leaves_cluster[:, 1] == i)]] for i in big_clusters]
#
# [len(i) for i in smiles_per_cluster]
# len(smiles_per_cluster)
#
#
# cluster_0 = smiles_per_cluster[0]
#
# D_0 = get_tanimoto_matrix(cluster_0, scaffolds=False, as_vector=True)
# tree_0 = hierarchy.complete(D_0)
#
# hierarchy.dendrogram(tree_0, count_sort=True, color_threshold=0.4, leaf_rotation=0)
# plt.show()
#
#
#
#
#
# # Full Tani = 0.8 & size = 100 -> 81 clusters
# # Full Tani = 0.4 & size = 100 -> 29 clusters
#
#
# # start out with a 0.8 pruned tani tree
# # prune that tree further with a 0.4 Tani threshold
#
#
# # start with a 0.8 pruned tree
# # take every big cluster (200>) and make a pruned tree out of it
# # makes the big tree clusters a superset of the small tree clusters
#
#
#
# smiles_per_cluster[2]
#
# clust_expanded  =[i for i in smiles_per_cluster if smi_2 in i]
# len(clust_expanded[0])
#
# smi_2 = 'Cc1cccc(Nc2nc(NCc3ccco3)nc(N)c2[N+](=O)[O-])c1'
#
# # Full Tani = 0.75 & size = 200 -> 21 clusters
# # Full Tani = 0.5 & size = 200 -> 11 clusters
# # Full Tani = 0.4 & size = 200 -> 5 clusters
#
# # scaf Tani = 0.25 & size = 200 -> 27 clusters
# # scaf Tani = 0.15 & size = 200 -> 24 clusters
# # scaf Tani = 0.5 & size = 200 -> 40 clusters
#
#
#
#
# D = torch.load('data/screen/tanimoto_matrix')
# D = torch.load('data/screen/tanimoto_scaffold_matrix')
#
# for i in tqdm(range(100000)):
#     smiles = smiles_all[np.random.choice(len(smiles_all), size=20, replace=False)]
#
#     D = get_tanimoto_matrix(smiles, verbose=False)
#     np.fill_diagonal(D, 0)
#     flat = np.sort(np.tril(D).flatten())
#
#     if np.min(flat[-5:]) > 0.6:
#         break
#
# from sklearn.cluster import SpectralClustering
# sc = SpectralClustering(50, affinity='precomputed', n_init=100)
# sc.fit_predict(S)
#
# sc.labels_
# Counter(sc.labels_)
#
# smiles_clust1 = smiles_all[np.where(sc.labels_ == 4)]
# clust_means = [np.mean(get_tanimoto_matrix(smiles_all[np.where(sc.labels_ == i)], verbose=False, scaffolds=True)) for i in range(50)]
#
#
# # df = pd.read_csv("tsne_perp50_500.csv")
# # df = df[df['split'] == 'test']
# # df['cluster'] = sc.labels_
# # df.to_csv("tsne_perp50_500_test.csv")
#
# clust = AgglomerativeClustering(metric="precomputed", linkage='complete', compute_full_tree=True, n_clusters=20,
#                                 compute_distances=True)
# clust.fit(1-D)
# linkage = linkage_matrix(clust)
# D_backup = D
#
#
# D = 1-D
# np.fill_diagonal(D, 0)
#
# linkage = hierarchy.complete(squareform(D))
#
# cut_clusters = hierarchy.cut_tree(linkage, n_clusters=[1, 3, 5])    # [sim_leaves]
# cut_clusters = hierarchy.cut_tree(linkage, height=[0.75])    # [sim_leaves]
#
# hierarchy.dendrogram(linkage, count_sort=True, color_threshold=0.7, labels=cut_clusters[:, -1], leaf_rotation=0)
# plt.show()
#
# cut = np.array((np.array(range(len(D))), cut_clusters[:, -1], )).T
# nodes_dist_threshold = linkage[linkage[:, 2] < 0.75][:, :2].flatten()
# leaves_dist_threshold = nodes_dist_threshold[nodes_dist_threshold < len(D)]
#
# leaves_cluster = cut[[True  if i in leaves_dist_threshold else False for i in cut[:, 0]]]
# # rebase the cluster labels
# for i, c in enumerate(set(leaves_cluster[:, 1])):
#     leaves_cluster[:, 1][leaves_cluster[:, 1] == c] = i
#
# big_clusters = [k for k, v in Counter(leaves_cluster[:, 1]).items() if v > 50]
#
#
#
#
# similarity_threshold = 0.4
#
# distance = 1 - similarity_threshold
# children = np.column_stack((np.array(range(len(clust.children_))).T, clust.children_))
# sim_children = children[clust.distances_ < distance]
# sim_leaves = sim_children[:, 1:][sim_children[:, 1:] < len(D)]
#
#
# cut_clusters = hierarchy.cut_tree(linkage, height=[0.9, 0.5, 0.1])    # [sim_leaves]
# clust.labels_.shape
#
#
# np.array((cut_clusters, np.array(range(len(smiles))))).T
#
# cut_children = sim_children[cut_nodes_idx]
# cut_leaves = sim_children[:, 1:][sim_children[:, 1:] < len(D)]
#
#
# cut_nodes_idx = hierarchy.cut_tree(linkage)
#
# cut_nodes_idx.shape
# D.shape
#
# cut_smiles = smiles_all[cut_nodes_idx]
# cut_clusters = clust.labels_[cut_nodes_idx]
#
# len(cut_smiles)
# len(smiles_all)
#
# possible_picks = smiles_all[sim_leaves]
#
#
#
# unique, counts = np.unique(clade_clusters, return_counts=True)
# big_clusters = unique[counts >= min_n_mols_per_clust]
# big_clusters_bool = [True if i in big_clusters else False for i in clade_clusters]
#
# clade_clusters = clade_clusters[big_clusters_bool]
# possible_picks = possible_picks[big_clusters_bool]
#
#
#
#
# get_tanimoto_matrix(possible_picks[clade_clusters == 0], verbose=False)
#
#
# unique, counts = np.unique(clade_clusters, return_counts=True)
#
#          # >= min_n_mols_per_clust
#
#
#
#
# # find the largest numbered node. Follow this node down and label every node along the way
# # repeat untill all nodes are labeled
#
# # get the top node of all nodes that have not been labeled
# # label the node with the current cluster number.
# # If the previous node was not all pure, go to the next connecting node
#
# # repeat untill al leaves are labeled
#
#
# node_clusers = np.zeros((len(sim_children)))
# cluster_count = 1
#
# while min(node_clusers) == 0:
#     top_node = sim_children[np.where((sim_children[:, 2] == np.max(sim_children[:, 2])) & (node_clusers == 0))[0]].flatten()
#     node_clusers[top_node[0]] = cluster_count
#
#     while top_node[2] < n:
#         # go to the next node
#         # set the cluster of this node
#         # update top node
#         connecting_node_idx = top_node[2] - n
#
#         node_clusers[sim_children[:, 0] == connecting_node_idx] = cluster_count
#
#         node_clusers[top_node[0]] = cluster_count
#         if top_node[2] >= n:
#             # check for leaves:
#             connecting_node_idx = top_node[2] - n
#             node_clusers[np.where(sim_children[:, 0] == connecting_node_idx)] = cluster_count
#
#
#
#
# sim_children[2][1]-20
#
# clust.distances_[clust.distances_ < 1 - similarity_threshold]
#
# # 11, 0, 18        5, 9       2, 7
# # sim_children
#
#
#
#
#
#
# D2 = get_tanimoto_matrix(possible_picks, verbose=False)
#
# D[6,9]
#
# clust.children_[5-len(D)]
# clust.distances_
#
#
#
#
# ii = itertools.count(D.shape[0])
# [{'node_id': next(ii), 'left': x[0], 'right':x[1]} for x in clust.children_]
#
# dict(enumerate(clust.children_, clust.n_leaves_))
#
# # The children of each non-leaf node. Values less than n_samples correspond to leaves of the tree which are the
# # original samples. A node i greater than or equal to n_samples is a non-leaf node and has
# # children children_[i - n_samples]. Alternatively at the i-th iteration, children[i][0] and children[i][1] are
# # merged to form node n_samples + i.
#
#
# from scipy.cluster.hierarchy import dendrogram
#
#
# def plot_dendrogram(clust, **kwargs):
#     # Create linkage matrix and then plot the dendrogram
#
#     # create the counts of samples under each node
#     counts = np.zeros(clust.children_.shape[0])
#     n_samples = len(clust.labels_)
#     for i, merge in enumerate(clust.children_):
#         current_count = 0
#         for child_idx in merge:
#             if child_idx < n_samples:
#                 current_count += 1  # leaf node
#             else:
#                 current_count += counts[child_idx - n_samples]
#         counts[i] = current_count
#
#     linkage_matrix = np.column_stack(
#         [clust.children_, clust.distances_, counts]
#     ).astype(float)
#
#     # Plot the corresponding dendrogram
#     dendrogram(linkage_matrix, truncate_mode="level", p=3)
#
# from matplotlib import pyplot as plt
# plot_dendrogram(clust, truncate_mode="level", p=3)
# plt.show()
#
#
#
#
#
# smiles = np.array(['Cc1c(NC(=O)COC(=O)Cc2ccc(Cl)cc2)c(=O)n(-c2ccccc2)n1C',
#        'Cc1c(C2C(C#N)=C([NH3+])N(c3cccc([N+](=O)[O-])c3)C3=C2C(=O)CCC3)cnn1C',
#        'CCC(=O)n1nc(NCc2ccc(OC)c(OC)c2)nc1NCc1ccc(OC)c(OC)c1',
#        'O=C(Nc1ccccc1)c1ccc2snnc2c1', 'COc1cccc(OCCSc2nc(N)cc(N)n2)c1',
#        'CCOc1ccc(NC(=O)c2cc3c(C)nc4ccccc4c3o2)cc1',
#        'COc1cc(NC(=O)Cc2nn(C)c(=O)c3ccccc23)cc(OC)c1',
#        'COc1ccc(Cl)cc1NC(=O)C(C)Sc1nncn1C',
#        'CC(C(=O)Nc1ccccc1C(=O)N1CCOCC1)N(c1ccccc1)S(C)(=O)=O',
#        'Cc1nc2ccccc2c2oc(C(=O)NCc3ccc4c(c3)OCO4)cc12',
#        'CC(C)CNC(=O)c1cc(-c2ccco2)on1',
#        'Cc1cc([N+](=O)[O-])nn1CCCC(=O)Nc1c(C)n(C)n(-c2ccccc2)c1=O',
#        'O=C(c1csnn1)N1CCC2(CCCN(c3ccc(-c4ccccc4)cc3)C2)CC1',
#        'Cc1ccc(CNC(=O)c2csc3c2CCCCC3)cc1',
#        'O=C1NCC2=C1C(c1ccccc1)Nc1ccccc1N2', 'CC(C)NC(=O)NC(=O)CSc1ncccn1',
#        'C(=N/N1CCN(c2ccccn2)CC1)\\c1ccco1',
#        'O=C(CSc1nc2ccccc2[nH]1)Nc1c(F)c(F)c(F)c(F)c1F',
#        'Cc1c(NC(=O)CN2C(=O)NC(CCc3ccccc3)C2=O)c(=O)n(-c2ccccc2)n1C',
#        'CSc1cccc(NC(=O)CN(c2cc(C)ccc2C)S(C)(=O)=O)c1'])