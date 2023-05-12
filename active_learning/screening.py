
import pandas as pd
import numpy as np
from active_learning.nn import Ensemble
from active_learning.data_prep import MasterDataset
from active_learning.data_handler import Handler
from active_learning.utils import Evaluate, to_torch_dataloader
from active_learning.acquisition import Acquisition, logits_to_pred
from tqdm.auto import tqdm
import torch
from math import ceil


def active_learning(n_start: int = 64, acquisition_method: str = 'exploration', max_screen_size: int = None,
                    batch_size: int = 16, architecture: str = 'gcn', seed: int = 0, bias: str = 'random',
                    optimize_hyperparameters: bool = False, ensemble_size: int = 5):

    representation = 'ecfp' if architecture == 'mlp' else 'graph'

    ds_screen = MasterDataset('screen', representation=representation)
    ds_test = MasterDataset('test', representation=representation)

    handler = Handler(n_start=n_start, seed=seed, bias=bias)
    ACQ = Acquisition(method=acquisition_method, seed=seed)
    eval_test = Evaluate()
    eval_screen = Evaluate()
    eval_train = Evaluate()
    x_test, y_test, smiles_test = ds_test.all()
    test_loader = to_torch_dataloader(x_test, y_test, batch_size=512, num_workers=4, shuffle=False, pin_memory=True)
    max_screen_size = len(ds_screen) if max_screen_size is None else max_screen_size

    hits_discovered, total_mols_screened, all_train_smiles = [], [], []

    for cycle in tqdm(range(ceil((max_screen_size - n_start)/batch_size)+1)):

        train_idx, screen_idx = handler()

        x_train, y_train, smiles_train = ds_screen[train_idx]
        x_screen, y_screen, smiles_screen = ds_screen[screen_idx]
        train_loader = to_torch_dataloader(x_train, y_train, batch_size=32, num_workers=4, shuffle=False, pin_memory=True)
        screen_loader = to_torch_dataloader(x_screen, y_screen, batch_size=512, num_workers=4, shuffle=False, pin_memory=True)

        all_train_smiles.append(';'.join(smiles_train.tolist()))
        hits_discovered.append(sum(y_train))
        hits = smiles_train[np.where(y_train == 1)]
        total_mols_screened.append(len(y_train))

        if len(train_idx) >= max_screen_size:
            break

        # we weigh classes based on the information in the training set.
        class_weights = torch.tensor([1 - sum((y_train == 0) * 1)/len(y_train),
                                      1 - sum((y_train == 1) * 1)/len(y_train)])

        M = Ensemble(seed=seed, ensemble_size=ensemble_size, architecture=architecture, class_weights=class_weights)

        if cycle == 0 and optimize_hyperparameters:
            M.optimize_hyperparameters(x_train, y_train, class_weights=class_weights)

        print("Training model")
        M.train(train_loader, verbose=False)

        print("Train/test/screen inference")
        train_logits_N_K_C = M.predict(train_loader)
        eval_train.eval(train_logits_N_K_C, y_train)

        test_logits_N_K_C = M.predict(test_loader)
        eval_test.eval(test_logits_N_K_C, y_test)

        screen_logits_N_K_C = M.predict(screen_loader)
        eval_screen.eval(screen_logits_N_K_C, y_screen)

        if len(train_idx) + batch_size > max_screen_size:
            batch_size = max_screen_size - len(train_idx)

        print("Sample acquisition")
        picks = ACQ.acquire(screen_logits_N_K_C, smiles_screen, hits=hits, n=batch_size)

        handler.add(picks)

    train_results = eval_train.to_dataframe("train_")
    test_results = eval_test.to_dataframe("test_")
    screen_results = eval_screen.to_dataframe('screen_')
    results = pd.concat([train_results, test_results, screen_results], axis=1)
    results['hits_discovered'] = hits_discovered
    results['total_mols_screened'] = total_mols_screened
    results['all_train_smiles'] = all_train_smiles

    return results
