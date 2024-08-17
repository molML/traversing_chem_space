
from active_learning.data_prep import smi_to_scaff
from tqdm.auto import tqdm
import torch
from scipy.spatial.distance import squareform
from rdkit.Chem import AllChem, DataStructs, Descriptors
import os
from active_learning.utils import eval_smarts
from warnings import warn
from active_learning.data_prep import MasterDataset
from rdkit import Chem
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
from rdkit.Chem import Descriptors
import umap


def pattern_finder(smiles: list[str], return_dict: bool = False):
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    patterns = torch.stack([mol_has_pattern(m, eval_smarts) for m in mols])

    if return_dict:
        return {k: v.item() for k, v in zip(eval_smarts.keys(), torch.sum(patterns, 0))}

    return torch.sum(patterns, 0)


def smiles_to_bitvec(smiles: str, scaffold: bool = False, radius: int = 2, nBits: int = 1024):
    if scaffold:
        return AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi_to_scaff(smiles, includeChirality=False)),
                                                     radius=radius, nBits=nBits)
    return AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), radius=radius, nBits=nBits)


def tani_matrix_w_db(smiles: list[str], bitvec_db: dict, as_vector: bool = True):

    smi_len = len(smiles)
    m = np.zeros([smi_len, smi_len], dtype=np.float16)  # We use 16-bit floats to prevent giant matrices
    # Calculate upper triangle of matrix
    for i in tqdm(range(smi_len), disable=True):
        for j in range(i, smi_len):
            m[i, j] = DataStructs.TanimotoSimilarity(bitvec_db[smiles[i]], bitvec_db[smiles[j]])
    # Fill in the lower triangle without having to loop (saves ~50% of time)
    m = m + m.T - np.diag(np.diag(m))
    # Fill the diagonal with 0's
    np.fill_diagonal(m, 0)
    if as_vector:
        return squareform(m)
    return m


def mols_to_descriptors(mols: list, progressbar: bool = False, normalize: bool = True) -> np.ndarray:
    """ Get the full set of available RDKit descriptors for a list of RDKit molecule objects

    :param mols: list of RDKit mol objects, e.g., as obtained through smiles_to_mols()
    :param progressbar: toggles progressbar (default = False)
    :param normalize: toggles min-max normalization
    :return: Numpy Array of all RDKit descriptors
    """
    x = np.array([list(Descriptors.CalcMolDescriptors(m).values()) for m in tqdm(mols, disable=not progressbar)])
    if normalize:
        x = max_normalization(x)
        if np.isnan(x).any():
            warn("There were some nan-values introduced by 0-columns. Replaced all nan-values with 0")
            x = np.nan_to_num(x, nan=0)

    return x


def compute_physchem(mols: list):

    X = []
    for m in mols:
        weight = Descriptors.ExactMolWt(m)
        logp = Descriptors.MolLogP(m)
        h_bond_donor = Descriptors.NumHDonors(m)
        h_bond_acceptors = Descriptors.NumHAcceptors(m)
        rotatable_bonds = Descriptors.NumRotatableBonds(m)
        atoms = Chem.rdchem.Mol.GetNumAtoms(m)
        heavy_atoms = Chem.rdchem.Mol.GetNumHeavyAtoms(m)
        molar_refractivity = Chem.Crippen.MolMR(m)
        topological_polar_surface_area = Chem.QED.properties(m).PSA
        formal_charge = Chem.rdmolops.GetFormalCharge(m)
        rings = Chem.rdMolDescriptors.CalcNumRings(m)

        X.append(np.array([weight, logp, h_bond_donor, h_bond_acceptors, rotatable_bonds, atoms, heavy_atoms,
                           molar_refractivity, topological_polar_surface_area, formal_charge, rings]))

    return np.array(X)


def find_start_end(seed=2, bias='large', architecture='mlp', batch_size=64, n_start=64, dataset='ALDH1',
                   acquisition_method='random', hits_only=False, retrain=True):

    df_subset = df[(df['seed'] == seed) &
                   (df['bias'] == bias) &
                   (df['architecture'] == architecture) &
                   (df['batch_size'] == batch_size) &
                   (df['n_start'] == n_start) &
                   (df['dataset'] == dataset) &
                   (df['retrain'] == retrain) &
                   (df['acquisition_method'] == acquisition_method)]

    start_smiles = df_subset[df_subset['train_cycle'] == 0]['all_train_smiles'].tolist()[0].split(';')
    end_smiles = df_subset[df_subset['train_cycle'] == 15]['all_train_smiles'].tolist()[0].split(';')

    smiles_idx = np.zeros(len(ds_screen))

    smiles_idx[[ds_screen.smiles_index[smi] for smi in end_smiles]] = 1
    smiles_idx[[ds_screen.smiles_index[smi] for smi in start_smiles]] = 2

    if hits_only:
        smiles_idx[ds_screen.y == 0] = 0

    return smiles_idx


def max_normalization(x: np.ndarray) -> np.ndarray:
    """ Perform max normalization on a matrix x / x.max(axis=0), just like
    sklearn.preprocessing.normalize(x, axis=0, norm='max')

    :param x: array to be normalized
    :return: normalized array
    """
    return x / x.max(axis=0)


def mol_has_pattern(mol, smarts):
    x = torch.zeros(len(smarts))
    for i, pattern in enumerate(smarts.values()):
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        if len(matches) > 0:
            x[i] = 1
    return x


def random_baseline(total_hits_in_library: int, batch_size, library_size: int = 100000,
                    total_screening_budget: int = 1000, start_size: int = 64) -> list:
    hits_found_in_theory = [0]
    n_mols_screened = 0
    while n_mols_screened < total_screening_budget:
        n = start_size if n_mols_screened == 0 else batch_size
        if n + n_mols_screened > 1000:
            n = 1000 - n_mols_screened

        if n_mols_screened == 0:
            hits = 1 + (n-1) * total_hits_in_library / library_size
        else:
            hits = n * total_hits_in_library / library_size
        hits_found_in_theory.append(hits + hits_found_in_theory[-1])

        n_mols_screened += n
        total_hits_in_library -= hits
        library_size -= n

    return hits_found_in_theory[1:]

if __name__ == '__main__':

    tot_dataset_hits = {'ALDH1': 4986, 'PKM2': 223, 'VDR': 239, 'IDH1': 11, 'ADRB2': 5, 'OPRK1': 9, 'KAT2A': 54, 'FEN1': 100, 'GBA': 55}
    baselines = {d: {s: {k: random_baseline(tot_dataset_hits[d], k, start_size=s) for k in [64, 32, 16]} for s in [2, 4, 8, 16, 32, 64, 128, 256, 512]} for d in ['ALDH1', 'PKM2', 'VDR', 'IDH1', 'ADRB2', 'OPRK1', 'KAT2A', 'FEN1', 'GBA']}


    files = [f for f in os.listdir('results') if f.endswith('results.csv')]
    all_dataframes = []
    for f in files:
        df = pd.read_csv(f'results/{f}')
        if not 'dataset' in df.columns:
            df['dataset'] = 'ALDH1'
        all_dataframes.append(df)

    df = pd.concat(all_dataframes)


    data = {'PKM2': {}, 'VDR': {}, 'ALDH1': {}, 'IDH1': {}, 'ADRB2': {}, 'OPRK1': {}, 'FEN1': {}, 'KAT2A': {}, 'GBA': {}}
    for dataset in data.keys():
        data[dataset]['index_smiles'] = torch.load(os.path.join('data', dataset, 'screen', 'index_smiles'))
        data[dataset]['smiles_index'] = torch.load(os.path.join('data', dataset, 'screen', 'smiles_index'))
        data[dataset]['all_y'] = torch.load(os.path.join('data', dataset, 'screen', 'y'))
        data[dataset]['all_x'] = torch.load(os.path.join('data', dataset, 'screen', 'x'))
        data[dataset]['all_smiles'] = torch.load(os.path.join('data', dataset, 'screen', 'smiles'))
        data[dataset]['bitvec'] = {smi: smiles_to_bitvec(smi) for smi in data[dataset]['all_smiles']}
        data[dataset]['scaff'] = {smi: Chem.MolToSmiles(Chem.MolFromSmiles(smi_to_scaff(smi, includeChirality=False))) for smi in data[dataset]['all_smiles']}


    train_smiles_per_run = df['all_train_smiles'].tolist()
    dataset_per_run = df['dataset'].tolist()
    df['hits_discovered'] = [int(i.replace('tensor(', '').replace(')', '')) if type(i) is str else i for i in df['hits_discovered'].tolist()]

    cum_hits = []
    for cycle, n_hits in zip(df['screen_cycle'], df['hits_discovered']):
        if cycle == 0:
            prev_hits = 0
        cum_hits.append(n_hits - prev_hits)
        prev_hits = n_hits
    df['cum_hits_discovered'] = cum_hits

    # precompute the occurence of all molecular patterns for each unique molecule
    patterns = {}
    for i, train_smiles in tqdm(enumerate(train_smiles_per_run)):
        train_smiles = np.array(train_smiles.split(';'))
        for smi in train_smiles:
            if smi not in patterns:
                patterns[smi] = pattern_finder([smi])


    # find which patterns occur more often in hits compared to non-hits
    hit_patterns = {ds: np.array(pattern_finder(data[ds]['all_smiles'][data[ds]['all_y'] == 1])) for ds in ['ALDH1']}
    nonhit_patterns = {ds: np.array(pattern_finder(data[ds]['all_smiles'])) for ds in [ 'ALDH1']}   #  [data[ds]['all_y'] == 0]

    df_patterns = pd.DataFrame({'pattern': eval_smarts.keys(), 'hit_patterns': hit_patterns['ALDH1'], 'nonhit_patterns': nonhit_patterns['ALDH1']})
    df_patterns['ratio'] = df_patterns['hit_patterns'] / df_patterns['nonhit_patterns']
    df_patterns = df_patterns.sort_values('ratio', ascending=False).reset_index()

    df_patterns.to_csv('figures/data/pattern_occurence_ALDH1.csv')

    torch.save(patterns, 'patterns.pt')
    torch.save(df, 'df.pt')


    # Precompute the similarity matrices for the end of each experiment so we can subset this matrix later.
    full_sim_matrices = {}
    full_sim_matrices_temp = {}
    for i, (train_smiles, dataset) in tqdm(enumerate(zip(train_smiles_per_run, dataset_per_run)), total=len(train_smiles_per_run)):
        train_smiles = np.array(train_smiles.split(';'))

        if len(train_smiles) <= 64:
            j = i

        if len(train_smiles) == 1000:
            full_sim_matrices_temp[j] = tani_matrix_w_db(train_smiles, data[dataset]['bitvec'], as_vector=False)

        if i % 10000 == 0:
            torch.save(full_sim_matrices_temp, f'full_sim_matrices/full_sim_matrices_{i}.pt')
            full_sim_matrices = full_sim_matrices | full_sim_matrices_temp
            full_sim_matrices_temp = {}
    # add the remaining stuff
    full_sim_matrices = full_sim_matrices | full_sim_matrices_temp
    torch.save(full_sim_matrices_temp, f'full_sim_matrices/full_sim_matrices_{i}.pt')

    # load everything again
    full_sim_matrices = {}
    for filename in os.listdir("full_sim_matrices"):
        if filename.startswith('full_sim_matrices_'):
            full_sim_matrices = full_sim_matrices | torch.load(f'full_sim_matrices/{filename}')
    patterns = torch.load('patterns.pt')
    df = torch.load('df.pt')


    train_smiles_per_run = df['all_train_smiles'].tolist()
    dataset_per_run = df['dataset'].tolist()
    batch_size_per_run = df['batch_size'].tolist()
    start_size_per_run = df['n_start'].tolist()

    df2 = {'mean_total_sims': [],
           'mean_total_hits_sims': [],
           'n_unique_scaffolds': [],
           'mean_tani_per_batch': [],
           'mean_tani_batch_to_start_batch': [],
           'mean_tani_all_mols_to_start_batch': [],
           'unique_patterns': [],
           **{k: [] for k in eval_smarts.keys()},
           'hit_unique_patterns': [],
           **{'hit_' + k: [] for k in eval_smarts.keys()}}


    for i, (train_smiles, dataset, batch_size, start_size) in tqdm(enumerate(zip(train_smiles_per_run, dataset_per_run, batch_size_per_run, start_size_per_run))):

        train_smiles = np.array(train_smiles.split(';'))
        size = len(train_smiles)

        train_scaffolds = [data[dataset]['scaff'][smi] for smi in train_smiles]
        df2['n_unique_scaffolds'].append(len(set(train_scaffolds)))

        if size <= 64:
            j = i
            batch_size = size
            batch_nr = -1
            batch_smiles = train_smiles
        else:
            batch_smiles = train_smiles[size - batch_size:size - batch_size + batch_nr * batch_size]
        len(batch_smiles)
        batch_nr += 1

        sim = full_sim_matrices[j]

        # all total mols uptill now sims
        all_current_sims = sim[:size, :size]
        df2['mean_total_sims'].append(np.mean(squareform(all_current_sims)))

        # Current batch sims
        if size == start_size:
            current_batch_sims = all_current_sims
        else:
            current_batch_sims = sim[size-batch_size:size-batch_size+batch_nr*batch_size, size-batch_size:size-batch_size+batch_nr*batch_size]
        df2['mean_tani_per_batch'].append(np.mean(squareform(current_batch_sims)))

        # sim to starting batch
        sim_all_mols_to_start = all_current_sims[:start_size, start_size:] if size != start_size else all_current_sims
        df2['mean_tani_all_mols_to_start_batch'].append(np.mean(sim_all_mols_to_start))    # has nans

        sim_batch_to_start = all_current_sims[:start_size, size-batch_size:size-batch_size+batch_nr*batch_size] if size != start_size else all_current_sims
        df2['mean_tani_batch_to_start_batch'].append(np.mean(sim_batch_to_start))    # has nans

        # find hits
        index_of_batch_smiles = np.array([data[dataset]['smiles_index'][smi] for smi in train_smiles])

        hit_smiles = train_smiles[data[dataset]['all_y'][index_of_batch_smiles] == 1]
        hit_indices = np.where(np.in1d(train_smiles, hit_smiles))[0]

        # hit sims
        df2['mean_total_hits_sims'].append(np.mean(squareform(np.array([all_current_sims[i, hit_indices] for i in hit_indices]))))

        # explored chemistry
        patt = torch.sum(torch.stack([patterns[smi] for smi in train_smiles]), 0)
        patt = {k: v.item() for k, v in zip(eval_smarts.keys(), patt)}
        for k, v in patt.items():
            df2[k].append(v)
        df2['unique_patterns'].append(sum(np.array(list(patt.values())) > 0))

        hit_patt = torch.sum(torch.stack([patterns[smi] for smi in hit_smiles]), 0)
        hit_patt = {k: v.item() for k, v in zip(eval_smarts.keys(), hit_patt)}
        for k, v in hit_patt.items():
            df2['hit_' + k].append(v)
        df2['hit_unique_patterns'].append(sum(np.array(list(hit_patt.values())) > 0))

    df2 = pd.DataFrame(df2)
    df = df.reset_index()
    df_combined = pd.concat([df, df2], axis=1)


    enrichment = []
    for i in tqdm(range(len(df_combined))):
        hits_found = df_combined["hits_discovered"][i]

        dataset = df_combined["dataset"][i]
        n_start = df_combined["n_start"][i]
        batch_size = df_combined["batch_size"][i]
        cycle = df_combined["train_cycle"][i]
        n_control_hits = baselines[dataset][n_start][batch_size][cycle]

        enrichment.append(hits_found/n_control_hits)
    df_combined['enrichment'] = enrichment

    # df_combined.to_csv('processed_results.csv', index=False)
    df_combined = df_combined.drop('all_train_smiles', axis=1)
    df_combined.to_csv('figures/data/processed_results.csv', index=False)


    # property ridge plot data
    df_ridge = df[df['dataset'] == 'ALDH1']
    df_ridge = df_ridge[df_ridge['batch_size'] == 64]
    df_ridge = df_ridge[df_ridge['n_start'] == 64]

    df_ridge['acquisition_method'][df_ridge['retrain'] == False] = 'exploitation no retrain'
    df_ridge = df_ridge.reset_index()
    df_ridge_long = {'train_cycle': [], 'smiles': [], 'hit': [], 'seed': [], 'architecture': [], 'dataset': [], 'bias': [],
                     'acquisition_method': [], 'LogP': [], 'MolWt': [], 'HBA': [], 'HBD': [], 'TPSA': [], }

    unique_df_ridge_smiles = list(set(sum([i.split(';') for i in df_ridge['all_train_smiles'].tolist()], [])))
    mol_db = {smi: Chem.MolFromSmiles(smi) for smi in unique_df_ridge_smiles}
    prop_sb = {smi: {'LogP': Descriptors.MolLogP(m),
                     'MolWt': Descriptors.MolWt(m),
                     'HBA': Descriptors.NumHAcceptors(m),
                     'HBD': Descriptors.NumHDonors(m),
                     'TPSA': Descriptors.TPSA(m)} for smi, m in mol_db.items()}

    for i in range(len(df_ridge)):
        tr_cycl = df_ridge['train_cycle'][i]
        seed = df_ridge['seed'][i]
        architecture = df_ridge['architecture'][i]
        dataset = df_ridge['dataset'][i]
        bias = df_ridge['bias'][i]
        acquisition_method = df_ridge['acquisition_method'][i]
        for smi in df_ridge['all_train_smiles'][i].split(';'):
            df_ridge_long['train_cycle'].append(tr_cycl)
            df_ridge_long['smiles'].append(smi)
            df_ridge_long['seed'].append(seed)
            df_ridge_long['architecture'].append(architecture)
            df_ridge_long['dataset'].append(dataset)
            df_ridge_long['bias'].append(bias)
            df_ridge_long['acquisition_method'].append(acquisition_method)
            df_ridge_long['hit'].append(data[dataset]['all_y'][data[dataset]['smiles_index'][smi]])
            df_ridge_long['LogP'].append(prop_sb[smi]['LogP'])
            df_ridge_long['MolWt'].append(prop_sb[smi]['MolWt'])
            df_ridge_long['HBA'].append(prop_sb[smi]['HBA'])
            df_ridge_long['HBD'].append(prop_sb[smi]['HBD'])
            df_ridge_long['TPSA'].append(prop_sb[smi]['TPSA'])

    df_ridge_long = pd.DataFrame(df_ridge_long)
    df_ridge_long.to_csv('figures/data/properties_ridge.csv', index=False)


    ### UMAP
    for ds in ['VDR', 'PKM2', 'ALDH1']:

        ds_screen = MasterDataset('screen', representation='ecfp', dataset=ds)
        mols = [Chem.MolFromSmiles(smi) for smi in ds_screen.smiles]
        props = mols_to_descriptors(mols, progressbar=True)

        U = umap.UMAP(n_neighbors=10, min_dist=0.25, spread=1)
        embedding = U.fit_transform(props)

        df_umap = pd.DataFrame({'UMAP1': embedding[:, 0],
                                'UMAP2': embedding[:, 1],
                                'y': ds_screen.y})

        df_umap.to_csv(f'figures/data/UMAP_{ds}.csv')

    df_umap['random'] = find_start_end(seed=3, bias='small', architecture='mlp', dataset='ALDH1', acquisition_method='random', hits_only=True)
    df_umap['bald'] = find_start_end(seed=3, bias='small', architecture='mlp', dataset='ALDH1', acquisition_method='bald', hits_only=True)
    df_umap['similarity'] = find_start_end(seed=3, bias='small', architecture='mlp', dataset='ALDH1', acquisition_method='similarity', hits_only=True)
    df_umap['exploitation'] = find_start_end(seed=3, bias='small', architecture='mlp', dataset='ALDH1', acquisition_method='exploitation', hits_only=True)
    df_umap['exploration'] = find_start_end(seed=3, bias='small', architecture='mlp', dataset='ALDH1', acquisition_method='exploration', hits_only=True)
    df_umap['exploitation_static'] = find_start_end(seed=3, bias='small', architecture='mlp', dataset='ALDH1', acquisition_method='exploitation', retrain=False, hits_only=True)

    df_umap.to_csv('figures/data/fig4_abc.csv')
