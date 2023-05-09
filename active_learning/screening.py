
import pandas as pd
import numpy as np
from active_learning.nn import Ensemble
from active_learning.data_prep import MasterDataset
from active_learning.data_handler import Handler
from active_learning.utils import Evaluate, to_torch_dataloader
from active_learning.acquisition import Acquisition, logits_to_pred
from tqdm import tqdm
import torch
from math import ceil


def active_learning(n_start: int = 100, acquisition_method: str = 'exploration', max_screen_size: int = None,
                    batch_size: int = 16, architecture: str = 'gcn', seed: int = 0,
                    optimize_hyperparameters: bool = False, bias: str = 'random'):

    representation = 'ecfp' if architecture == 'mlp' else 'graph'

    ds_screen = MasterDataset('screen', representation=representation)
    ds_test = MasterDataset('test', representation=representation)

    handler = Handler(n_start=n_start, seed=seed, bias=bias)
    ACQ = Acquisition(method=acquisition_method, seed=seed)
    eval_test = Evaluate()
    eval_screen = Evaluate()
    x_test, y_test, smiles_test = ds_test.all()
    test_loader = to_torch_dataloader(x_test, y_test, batch_size=256, num_workers=4, shuffle=False)
    max_screen_size = len(ds_screen) if max_screen_size is None else max_screen_size

    hits_discovered, total_mols_screened, all_train_smiles = [], [], []

    for cycle in tqdm(range(ceil((max_screen_size - n_start)/batch_size)+1)):
        # break
        train_idx, screen_idx = handler()

        x_train, y_train, smiles_train = ds_screen[train_idx]
        x_screen, y_screen, smiles_screen = ds_screen[screen_idx]
        train_loader = to_torch_dataloader(x_train, y_train, batch_size=256, num_workers=4, shuffle=False)
        screen_loader = to_torch_dataloader(x_screen, y_screen, batch_size=256, num_workers=4, shuffle=False)

        all_train_smiles.append(';'.join(smiles_train.tolist()))
        hits_discovered.append(sum(y_train))
        hits = smiles_train[np.where(y_train == 1)]
        total_mols_screened.append(len(y_train))

        if len(train_idx) >= max_screen_size:
            break

        # we weigh classes based on the information in the training set.
        class_weights = torch.tensor([1 - sum((y_train == 0) * 1)/len(y_train),
                                      1 - sum((y_train == 1) * 1)/len(y_train)])

        M = Ensemble(seed=seed, ensemble_size=3, architecture=architecture, epochs=100, l2_lambda=1e-6,
                     class_weights=class_weights)

        if cycle == 0 and optimize_hyperparameters:
            M.optimize_hyperparameters(x_train, y_train)

        print("Training model")
        M.train(train_loader, verbose=False)

        print("Test inference")
        test_logits_N_K_C = M.predict(test_loader)
        y_hat_test = logits_to_pred(test_logits_N_K_C, return_uncertainty=False, return_prob=False)
        eval_test.eval(y_hat_test, y_test)

        print("Screen inference")
        screen_logits_N_K_C = M.predict(screen_loader)
        y_hat_screen = logits_to_pred(screen_logits_N_K_C, return_uncertainty=False, return_prob=False)
        eval_screen.eval(y_hat_screen, y_screen)

        if len(train_idx) + batch_size > max_screen_size:
            batch_size = max_screen_size - len(train_idx)

        print("Sample acquisition")
        picks = ACQ.acquire(screen_logits_N_K_C, smiles_screen, hits=hits, n=batch_size)

        handler.add(picks)

    test_results = eval_test.to_dataframe("test_")
    screen_results = eval_screen.to_dataframe('screen_')
    results = pd.concat([test_results, screen_results], axis=1)
    results['hits_discovered'] = hits_discovered
    results['total_mols_screened'] = total_mols_screened
    results['all_train_smiles'] = all_train_smiles

    return results
