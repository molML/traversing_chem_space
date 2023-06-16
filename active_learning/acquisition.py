"""

This script contains a selection of sample acquisition methods for active learning.
All functions here operate on model predictions in the form of logits_N_K_C = [N, num_inference_samples, num_classes].
Here N is a molecule, K are the number of sampled predictions (i.e., 10 for a 10-model ensemble), and C = 2 ([0, 1]):

    - Acquisition: class that handles all molecule acquisition
    - logits_to_pred: Get the probabilities/class vector and sample uncertainty from logits
    - logits_mean: Get the logit mean with the logsumexp trick
    - entropy: Calculates the Shannon Entropy on the logits
    - mean_sample_entropy: Calculates the mean entropy for each sample given multiple ensemble predictions
    - mutual_information: Calculates the Mutual Information
    - greedy_exploitation: Get the n highest predicted samples
    - greedy_exploration: Get the n most uncertain samples (based on entropy)
    - dynamic_exploration: Gradually move from greedy exploration to greedy exploitation
    - bald: Get the n molecules with the lowest Mutual Information - Houlsby et al., 2011
    - batch_bald: Get BatchBALD batch - Kirch et al., 2019, NeurIPS
    - similarity_search: Perform similarity search, take the n screen mols with the highest similarity to any hit

    Author: Derek van Tilborg, Eindhoven University of Technology, May 2023

"""

import numpy as np
import torch
from torch import Tensor, tensor
from batchbald_redux.batchbald import get_batchbald_batch
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect as ECFPbitVec
from rdkit.DataStructs import BulkTanimotoSimilarity
from rdkit import Chem
import math


class Acquisition:
    def __init__(self, method: str, seed: int = 42, **kwargs):

        self.acquisition_method = {'random': self.random_pick,
                                   'exploration': greedy_exploration,
                                   'exploitation': greedy_exploitation,
                                   'dynamic': dynamic_exploration,
                                   'dynamic_bald': dynamic_exploration_bald,
                                   'bald': bald,
                                   'batch_bald': batch_bald,
                                   'similarity': similarity_search}

        assert method in self.acquisition_method.keys(), f"Specified 'method' not available. " \
                                                         f"Select from: {self.acquisition_method.keys()}"

        self.method = method
        self.params = kwargs
        self.rng = np.random.default_rng(seed=seed)
        self.iteration = 0

    def acquire(self, logits_N_K_C: Tensor, smiles: np.ndarray[str], hits: np.ndarray[str], n: int = 1) -> \
            np.ndarray[str]:

        self.iteration += 1

        return self.acquisition_method[self.method](logits_N_K_C=logits_N_K_C, smiles=smiles, n=n, hits=hits,
                                                    iteration=self.iteration, **self.params)

    def __call__(self, *args, **kwargs) -> np.ndarray[str]:
        return self.acquire(*args, **kwargs)

    def random_pick(self, smiles: np.ndarray[str], n: int = 1, return_smiles: bool = True, **kwargs) -> np.ndarray:
        """ select n random samples """
        picks_idx = self.rng.integers(0, len(smiles), n)

        return smiles[picks_idx] if return_smiles else picks_idx


def logits_to_pred(logits_N_K_C: Tensor, return_prob: bool = True, return_uncertainty: bool = True) -> (Tensor, Tensor):
    """ Get the probabilities/class vector and sample uncertainty from the logits """

    mean_probs_N_C = torch.mean(torch.exp(logits_N_K_C), dim=1)
    uncertainty = mean_sample_entropy(logits_N_K_C)

    if return_prob:
        y_hat = mean_probs_N_C
    else:
        y_hat = torch.argmax(mean_probs_N_C, dim=1)

    if return_uncertainty:
        return y_hat, uncertainty
    else:
        return y_hat


def logit_mean(logits_N_K_C: Tensor, dim: int, keepdim: bool = False) -> Tensor:
    """ Logit mean with the logsumexp trick - Kirch et al., 2019, NeurIPS """

    return torch.logsumexp(logits_N_K_C, dim=dim, keepdim=keepdim) - math.log(logits_N_K_C.shape[dim])


def entropy(logits_N_K_C: Tensor, dim: int, keepdim: bool = False) -> Tensor:
    """Calculates the Shannon Entropy """

    return -torch.sum((torch.exp(logits_N_K_C) * logits_N_K_C).double(), dim=dim, keepdim=keepdim)


def mean_sample_entropy(logits_N_K_C: Tensor, dim: int = -1, keepdim: bool = False) -> Tensor:
    """Calculates the mean entropy for each sample given multiple ensemble predictions - Kirch et al., 2019, NeurIPS"""

    sample_entropies_N_K = entropy(logits_N_K_C, dim=dim, keepdim=keepdim)
    entropy_mean_N = torch.mean(sample_entropies_N_K, dim=1)

    return entropy_mean_N


def mutual_information(logits_N_K_C: Tensor) -> Tensor:
    """ Calculates the Mutual Information - Kirch et al., 2019, NeurIPS """

    # this term represents the entropy of the model prediction (high when uncertain)
    entropy_mean_N = mean_sample_entropy(logits_N_K_C)

    # This term is the expectation of the entropy of the model prediction for each draw of model parameters
    mean_entropy_N = entropy(logit_mean(logits_N_K_C, dim=1), dim=-1)

    I = mean_entropy_N - entropy_mean_N

    return I


def greedy_exploitation(logits_N_K_C: Tensor, smiles: np.ndarray[str], n: int = 1, **kwargs) -> np.ndarray[str]:
    """ Get the n highest predicted samples """

    mean_probs_hits = torch.mean(torch.exp(logits_N_K_C), dim=1)[:, 1]
    picks_idx = torch.argsort(mean_probs_hits, descending=True)[:n]

    return np.array([smiles[picks_idx.cpu()]]) if n == 1 else smiles[picks_idx.cpu()]


def greedy_exploration(logits_N_K_C: Tensor, smiles: np.ndarray[str], n: int = 1, **kwargs) -> np.ndarray[str]:
    """ Get the n most samples with the most variance in hit classification """

    entropy_mean_N = mean_sample_entropy(logits_N_K_C)
    # sd_mean_N = torch.std(torch.exp(logits_N_K_C), dim=1)[:, 1]
    picks_idx = torch.argsort(entropy_mean_N, descending=True)[:n]

    return np.array([smiles[picks_idx.cpu()]]) if n == 1 else smiles[picks_idx.cpu()]


def dynamic_exploration(logits_N_K_C: Tensor, smiles: np.ndarray[str], n: int = 1, lambd: float = 0.95,
                        iteration: int = 0, **kwargs) -> np.ndarray[str]:
    """ starts with 100% exploration, approaches the limit of 100% exploitation. The speed in which we stop
    exploring depends on lambda. For example, a lambda of 0.9 will require 44 iterations to reach full exploitation,
    a lambda of 0.5 will get there in only 7 iterations """

    exploitation_factor = (1/(lambd ** iteration)) - 1
    n_exploit = round(n * exploitation_factor)
    n_explore = n - n_exploit

    exploitative_picks = greedy_exploitation(logits_N_K_C, smiles, n=n_exploit)
    explorative_picks = greedy_exploration(logits_N_K_C, smiles, n=n_explore)

    return np.concatenate((exploitative_picks, explorative_picks))


def dynamic_exploration_bald(logits_N_K_C: Tensor, smiles: np.ndarray[str], n: int = 1, lambd: float = 0.95,
                             iteration: int = 0, **kwargs) -> np.ndarray[str]:
    """ starts with 100% exploration, approaches the limit of 100% exploitation. The speed in which we stop
    exploring depends on lambda. For example, a lambda of 0.9 will require 44 iterations to reach full exploitation,
    a lambda of 0.5 will get there in only 7 iterations """

    exploitation_factor = (1/(lambd ** iteration)) - 1
    n_exploit = round(n * exploitation_factor)
    n_explore = n - n_exploit

    exploitative_picks = greedy_exploitation(logits_N_K_C, smiles, n=n_exploit)
    explorative_picks = bald(logits_N_K_C, smiles, n=n_explore)

    return np.concatenate((exploitative_picks, explorative_picks))


def bald(logits_N_K_C: Tensor, smiles: np.ndarray[str], n: int = 1, **kwargs) -> np.ndarray[str]:
    """ Get the n molecules with the lowest Mutual Information """
    I = mutual_information(logits_N_K_C)

    picks_idx = torch.argsort(I, descending=False)[:n]

    return smiles[picks_idx.cpu()]


def batch_bald(logits_N_K_C: Tensor, smiles: np.ndarray[str], n: int = 1, num_samples: int = 1000, **kwargs) -> \
        np.ndarray[str]:
    """ Get BatchBALD batch - Kirch et al., 2019, NeurIPS"""
    if n == 1:
        return bald(logits_N_K_C, smiles)

    num_samples = logits_N_K_C.shape[0] if num_samples is None else num_samples
    candidate_batch = get_batchbald_batch(logits_N_K_C, batch_size=n, num_samples=num_samples,
                                          device=torch.device("cuda" if torch.cuda.is_available() else "cpu"))

    return smiles[candidate_batch.indices]


def similarity_search(hits: np.ndarray[str], smiles: np.ndarray[str], n: int = 1, radius: int = 2, nBits: int = 1024,
                      **kwargs) -> np.ndarray[str]:
    """ 1. Compute the similarity of all screen smiles to all hit smiles
        2. take the n screen smiles with the highest similarity to any hit """

    fp_hits = [ECFPbitVec(Chem.MolFromSmiles(smi), radius=radius, nBits=nBits) for smi in hits]
    fp_smiles = [ECFPbitVec(Chem.MolFromSmiles(smi), radius=radius, nBits=nBits) for smi in smiles]

    m = np.zeros([len(hits), len(smiles)], dtype=np.float16)
    for i in range(len(hits)):
        m[i] = BulkTanimotoSimilarity(fp_hits[i], fp_smiles)

    # get the n highest similarity smiles to any hit
    picks_idx = np.argsort(np.max(m, axis=0))[::-1][:n]

    return smiles[picks_idx]
