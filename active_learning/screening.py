
import pandas as pd
from active_learning.nn.models import GCNEnsemble, BayesianNN
from active_learning.data_prep import MasterDataset
from active_learning.data_handler import Handler
from active_learning.utils import Evaluate
from active_learning.acquisition import Acquisition
from tqdm import tqdm
from math import ceil


def active_learning(n_start: int = 100, balanced_start: float = False, acquisition_method: str = 'greedy_explorative',
                    max_screen_size: int = None, batch_size: int = 10, architecture: str = 'gcn', seed: int = 0):

    representation = 'ecfp' if architecture == 'bnn' else 'graph'

    ds_screen = MasterDataset('screen', representation=representation)
    ds_test = MasterDataset('test', representation=representation)

    handler = Handler(n_train=n_start, balanced=balanced_start)
    ACQ = Acquisition(method=acquisition_method, seed=seed)
    eval_test = Evaluate()
    eval_screen = Evaluate()
    x_test, y_test, smiles_test = ds_test.all()
    max_screen_size = len(ds_screen) if max_screen_size is None else max_screen_size

    hits_discovered, total_mols_screened = [], []

    for cycle in tqdm(range(ceil((max_screen_size - n_start)/batch_size)+1)):

        train_idx, screen_idx = handler()

        x_train, y_train, smiles_train = ds_screen[train_idx]
        x_screen, y_screen, smiles_screen = ds_screen[screen_idx]

        hits_discovered.append(sum(y_train))
        total_mols_screened.append(len(y_train))

        if len(train_idx) >= max_screen_size:
            break

        if architecture == 'gcn':
            M = GCNEnsemble(seed=seed, ensemble_size=3, epochs=10)
        elif architecture == 'bnn':
            M = BayesianNN(seed=seed, to_gpu=False, epochs=100)

        # if cycle == 0:
        #     M.optimize_hyperparameters(x_train, y_train)

        M.train(x_train, y_train)

        y_hat, y_hat_mu, y_hat_sigma = M.predict(x_test)
        eval_test.eval(y_hat_mu, y_test)

        y_hat, y_hat_mu, y_hat_sigma = M.predict(x_screen)
        eval_screen.eval(y_hat_mu, y_screen)

        if len(train_idx) + batch_size > max_screen_size:
            batch_size = max_screen_size - len(train_idx)

        picks = ACQ.acquire(y_hat_mu, y_hat_sigma, smiles_screen, n=batch_size)

        handler.add(picks)

    test_results = eval_test.to_dataframe("test_")
    screen_results = eval_screen.to_dataframe('screen_')
    results = pd.concat([test_results, screen_results], axis=1)
    results['hits_discovered'] = hits_discovered
    results['total_mols_screened'] = total_mols_screened

    return results