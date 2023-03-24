
import numpy as np
from active_learning.utils import smiles_to_ecfp
from rdkit.DataStructs import BulkTanimotoSimilarity
import torch
ACQUISITION_METHODS = ['random', 'greedy_exploitative', 'greedy_explorative', 'epsilon_greedy_exploitative',
                          'epsilon_greedy_explorative', 'diversity_explorative', 'diversity_exploitative',
                          'dynamic_uncertainty', 'dynamic_diversity']


class Acquisition:
    def __init__(self, method: str, seed: int = 42, **kwargs):
        assert method in ACQUISITION_METHODS, f"Specified 'method' not available. Select from: {ACQUISITION_METHODS}"

        self.method = method
        self.params = kwargs
        self.rng = np.random.default_rng(seed=seed)
        self.iteration = 0

    def acquire(self, y_hat_mu: np.ndarray, y_hat_sigma: np.ndarray, smiles: np.ndarray, n: int = 1) -> np.ndarray:

        self.iteration += 1

        if self.method == 'random':
            return self.random_pick(smiles=smiles, n=n)
        if self.method == 'greedy_exploitative':
            return self.greedy_exploitative(y_hat_mu=y_hat_mu, smiles=smiles, n=n)
        if self.method == 'greedy_explorative':
            return self.greedy_explorative(y_hat_sigma=y_hat_sigma, smiles=smiles, n=n)
        if self.method == 'epsilon_greedy_exploitative':
            return self.epsilon_greedy_exploitative(y_hat_mu=y_hat_mu, smiles=smiles, n=n, **self.params)
        if self.method == 'epsilon_greedy_explorative':
            return self.epsilon_greedy_explorative(y_hat_sigma=y_hat_sigma, smiles=smiles, n=n, **self.params)
        if self.method == 'diversity_exploitative':
            return self.diversity_exploitative(y_hat_mu=y_hat_mu, smiles=smiles, n=n)
        if self.method == 'diversity_explorative':
            return self.diversity_explorative(y_hat_sigma=y_hat_sigma, smiles=smiles, n=n)
        if self.method == 'dynamic_uncertainty':
            return self.dynamic_uncertainty(y_hat_mu=y_hat_mu, y_hat_sigma=y_hat_sigma, smiles=smiles, n=n,
                                            **self.params)
        if self.method == 'dynamic_diversity':
            return self.dynamic_diversity(y_hat_mu=y_hat_mu, y_hat_sigma=y_hat_sigma, smiles=smiles, n=n,
                                          **self.params)

    def __call__(self, y_hat_mu: np.ndarray, y_hat_sigma: np.ndarray, smiles: np.ndarray, n: int = 1, **kwargs) -> \
            np.ndarray:
        return self.acquire(y_hat_mu, y_hat_sigma, smiles, n, **kwargs)

    def random_pick(self, smiles: np.ndarray, n: int = 1, return_smiles: bool = True) -> np.ndarray:
        """ select n random samples """
        picks_idx = self.rng.integers(0, len(smiles), n)

        return smiles[picks_idx] if return_smiles else picks_idx

    @ staticmethod
    def greedy_exploitative(y_hat_mu: np.ndarray, smiles: np.ndarray, n: int = 1,
                            return_smiles: bool = True) -> np.ndarray:
        """ Get the n highest predicted samples """
        y_hat_mu = np.array(y_hat_mu.cpu()) if type(y_hat_mu) is torch.Tensor else y_hat_mu
        picks_idx = np.argsort(y_hat_mu)[::-1][:n]

        return smiles[picks_idx] if return_smiles else picks_idx

    @ staticmethod
    def greedy_explorative(y_hat_sigma: np.ndarray, smiles: np.ndarray, n: int = 1,
                           return_smiles: bool = True) -> np.ndarray:
        """ Get the n most uncertain samples """
        y_hat_sigma = np.array(y_hat_sigma.cpu()) if type(y_hat_sigma) is torch.Tensor else y_hat_sigma
        picks_idx = np.argsort(y_hat_sigma)[::-1][:n]

        return smiles[picks_idx] if return_smiles else picks_idx

    def epsilon_greedy_exploitative(self, y_hat_mu: np.ndarray, smiles: np.ndarray, n: int = 1,
                                    epsilon: float = 0.1) -> np.ndarray:
        """ greedy exploitative with an epsilon chance of selecting a random sample """

        n_random_picks = sum((self.rng.random(n) < epsilon) * 1)
        n_greedy_picks = n - n_random_picks

        greedy_picks = self.greedy_exploitative(y_hat_mu, smiles, n=n_greedy_picks)
        random_picks = self.random_pick(np.array([i for i in smiles if i not in greedy_picks]), n=n_random_picks)

        return np.concatenate((greedy_picks, random_picks))

    def epsilon_greedy_explorative(self, y_hat_sigma: np.ndarray, smiles: np.ndarray, n: int = 1,
                                   epsilon: float = 0.1) -> np.ndarray:
        """ greedy exploitative with an epsilon chance of selecting a random sample """

        n_random_picks = sum((self.rng.random(n) < epsilon) * 1)
        n_greedy_picks = n - n_random_picks

        greedy_picks = self.greedy_explorative(y_hat_sigma, smiles, n=n_greedy_picks)
        random_picks = self.random_pick(np.array([i for i in smiles if i not in greedy_picks]), n=n_random_picks)

        return np.concatenate((greedy_picks, random_picks))

    def clusters(self, y_hat: np.ndarray, y_hat_sigma: np.ndarray, n: int = 1):
        """ Select the n samples most representative of n clusters in the data """
        raise NotImplementedError

    def diversity_explorative(self, y_hat_sigma: np.ndarray, smiles: np.ndarray, n: int = 1,
                              return_smiles: bool = True) -> np.ndarray:
        """ select the n molecules that are structurally most similar to the n best predicted values """

        picks = self.greedy_explorative(y_hat_sigma, smiles, return_smiles=False, n=n)
        fps = smiles_to_ecfp(smiles, to_array=False)

        S = np.array([BulkTanimotoSimilarity(fps[pick], fps) for pick in picks])
        most_similar = np.argsort(S)[:, -2]  # we index -2 because -1 is the pick itself

        return smiles[most_similar] if return_smiles else most_similar

    def diversity_exploitative(self, y_hat_mu: np.ndarray, smiles: np.ndarray, n: int = 1,
                               return_smiles: bool = True) -> np.ndarray:
        """ select the n molecules that are structurally most similar to the n best predicted values """

        picks = self.greedy_exploitative(y_hat_mu, smiles, return_smiles=False, n=n)
        fps = smiles_to_ecfp(smiles, to_array=False)

        S = np.array([BulkTanimotoSimilarity(fps[pick], fps) for pick in picks])
        most_similar = np.argsort(S)[:, -2]  # we index -2 because -1 is the pick itself

        return smiles[most_similar] if return_smiles else most_similar

    def dynamic_uncertainty(self, y_hat_mu: np.ndarray, y_hat_sigma: np.ndarray, smiles: np.ndarray, n: int = 1,
                            lambd: float = 0.9) -> np.ndarray:
        """ starts with 100% exploration, approaches the limit of 100% exploitation. The speed in which we stop
        exploring depends on lambda. For example, a lambda of 0.9 will require 44 iterations to reach full exploitation,
        a lambda of 0.5 will get there in only 7 iterations """

        exploration_factor = ((lambd ** self.iteration) / 1)
        n_explore = round(n * exploration_factor)
        n_exploit = n - n_explore

        exploitative_picks = self.greedy_exploitative(y_hat_mu, smiles, n=n_exploit)
        explorative_picks = self.greedy_explorative(y_hat_sigma, smiles, n=n_explore)

        return np.concatenate((exploitative_picks, explorative_picks))

    def dynamic_diversity(self, y_hat_mu: np.ndarray, y_hat_sigma: np.ndarray, smiles: np.ndarray, n: int = 1,
                          lambd: float = 0.9) -> np.ndarray:

        exploration_factor = ((lambd ** self.iteration) / 1)
        n_explore = round(n * exploration_factor)
        n_exploit = n - n_explore

        exploitative_picks = self.diversity_exploitative(y_hat_mu, smiles, n=n_exploit)
        explorative_picks = self.diversity_explorative(y_hat_sigma, smiles, n=n_explore)

        return np.concatenate((exploitative_picks, explorative_picks))
