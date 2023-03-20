
import numpy as np
import torch


class Acquisition:
    def __init__(self, method: str, seed: int = 42, **kwargs):
        assert method in ['random', 'greedy_exploitative', 'greedy_explorative', 'epsilon_greedy_exploitative',
                          'epsilon_greedy_explorative']

        self.method = method
        self.params = kwargs
        self.rng = np.random.default_rng(seed=seed)
        self.iteration = 0

    def acquire(self, y_hat, y_hat_sigma, smiles: np.ndarray, n: int = 1, **kwargs):
        self.iteration += 1
        if self.method == 'random':
            return self.random_pick(smiles=smiles, n=n)
        if self.method == 'greedy_exploitative':
            return self.greedy_exploitative(y_hat=y_hat, smiles=smiles, n=n)
        if self.method == 'greedy_explorative':
            return self.greedy_explorative(y_hat_sigma=y_hat_sigma, smiles=smiles, n=n)
        if self.method == 'epsilon_greedy_exploitative':
            return self.epsilon_greedy_exploitative(y_hat=y_hat, smiles=smiles, n=n, **kwargs)
        if self.method == 'epsilon_greedy_explorative':
            return self.epsilon_greedy_explorative(y_hat_sigma=y_hat_sigma, smiles=smiles, n=n, **kwargs)

    def __call__(self, n, y_hat, y_hat_sigma, **kwargs):
        return self.acquire(n, y_hat, y_hat_sigma, **kwargs)

    def random_pick(self, smiles: np.ndarray, n: int = 1):
        """ select n random samples """
        picks_idx = self.rng.integers(0, len(smiles), n)

        return smiles[picks_idx]

    @ staticmethod
    def greedy_exploitative(y_hat, smiles: np.ndarray, n: int = 1):
        """ Get the n highest predicted samples """
        y_hat = np.array(y_hat) if type(y_hat) is torch.Tensor else y_hat
        picks_idx = np.argsort(y_hat)[::-1][:n]

        return smiles[picks_idx]

    @ staticmethod
    def greedy_explorative(y_hat_sigma, smiles: np.ndarray, n: int = 1):
        """ Get the n most uncertain samples """
        y_hat_sigma = np.array(y_hat_sigma) if type(y_hat_sigma) is torch.Tensor else y_hat_sigma
        picks_idx = np.argsort(y_hat_sigma)[::-1][:n]

        return smiles[picks_idx]

    def epsilon_greedy_exploitative(self, y_hat, smiles: np.ndarray, n: int = 1, epsilon: float = 0.1):
        """ greedy exploitative with an epsilon chance of selecting a random sample """

        n_random_picks = sum((self.rng.random(n) < epsilon) * 1)
        n_greedy_picks = n - n_random_picks

        greedy_picks = self.greedy_exploitative(y_hat, smiles, n=n_greedy_picks)
        random_picks = self.random_pick(np.array([i for i in smiles if i not in greedy_picks]), n=n_random_picks)

        return np.concatenate((greedy_picks, random_picks))

    def epsilon_greedy_explorative(self, y_hat_sigma, smiles: np.ndarray, n: int = 1, epsilon: float = 0.1):
        """ greedy exploitative with an epsilon chance of selecting a random sample """

        n_random_picks = sum((self.rng.random(n) < epsilon) * 1)
        n_greedy_picks = n - n_random_picks

        greedy_picks = self.greedy_explorative(y_hat_sigma, smiles, n=n_greedy_picks)
        random_picks = self.random_pick(np.array([i for i in smiles if i not in greedy_picks]), n=n_random_picks)

        return np.concatenate((greedy_picks, random_picks))

    def clusters(self, y_hat, y_hat_sigma, n: int = 1):
        """ Select the n samples most representative of n clusters in the data """
        pass
