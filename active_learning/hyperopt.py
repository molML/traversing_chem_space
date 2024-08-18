"""

Code for hyperparamter optimization using bootstrapped n-fold cross-validation.
    - optimize_hyperparameters()
    - k_fold_cross_validation()

Derek van Tilborg | 06-03-2023 | Eindhoven University of Technology

"""

import itertools
from operator import itemgetter
import numpy as np
import torch
from tqdm.auto import tqdm
from sklearn.metrics import balanced_accuracy_score
from active_learning.hyperparameters import MLP_hypers, GCN_hypers
from active_learning.utils import to_torch_dataloader
from active_learning.acquisition import logits_to_pred


def optimize_hyperparameters(x: np.ndarray, y: np.ndarray, class_weights=None, n_folds: int = 5,
                             architecture="gcn") -> dict:
    """ Function to optimize hyperparameters on a dataset using bootstrapped k-fold cross-validation """

    assert architecture in ['mlp', 'gcn'], f"'architecture' must be 'mlp' or 'gcn'"

    if architecture == 'gcn':
        hypers = GCN_hypers
    else:
        hypers = MLP_hypers

    all_hypers = [dict(zip(hypers.keys(), v)) for v in itertools.product(*hypers.values())]
    score_hypers = []
    for hyper in tqdm(all_hypers):
        hyper['ensemble_size'] = 1
        score = k_fold_cross_validation(x, y, n_folds=n_folds, seed=42, architecture=architecture, verbose=True,
                                        class_weights=class_weights, **hyper)
        score_hypers.append((score, hyper))

    best_hypers = sorted(score_hypers, key=itemgetter(0))[::-1][0][1]

    return best_hypers


def k_fold_cross_validation(x, y, n_folds: int = 5, seed: int = 42, architecture: str = 'mlp', verbose: bool = False,
                            class_weights=None, **kwargs) -> float:
    from active_learning.nn import Ensemble

    # Define some variables
    # Set random state and create folds
    rng = np.random.default_rng(seed)
    folds = rng.integers(low=0, high=n_folds, size=len(x))
    scores = []

    for i in range(n_folds):
        # Subset train/test folds
        if type(x) is np.ndarray:
            x_train, y_train = x[folds != i], y[folds != i]
            x_test, y_test = x[folds == i], y[folds == i]
        else:
            x_train, y_train = [x[j] for j in np.where(folds != i)[0]], y[folds != i]
            x_test, y_test = [x[j] for j in np.where(folds == i)[0]], y[folds == i]

        train_loader = to_torch_dataloader(x_train, y_train, batch_size=256, num_workers=4, shuffle=False)
        test_loader = to_torch_dataloader(x_test, y_test, batch_size=256, num_workers=4, shuffle=False)

        m = Ensemble(architecture=architecture, class_weights=class_weights, **kwargs)
        m.train(train_loader, verbose=verbose)

        test_logits_N_K_C = m.predict(test_loader)
        y_hat_test = logits_to_pred(test_logits_N_K_C, return_uncertainty=False, return_prob=False)
        scores.append(balanced_accuracy_score(y_test, y_hat_test.cpu()))

        # Delete the NN to free memory
        try:
            del m
            torch.cuda.empty_cache()
        except:
            pass

    return sum(scores) / len(scores)
