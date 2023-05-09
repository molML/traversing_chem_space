"""

Code for Bayesian hyperparamter optimization using bootstrapped n-fold cross-validation.
    - optimize_hyperparameters()
    - BayesianOptimization
    - dict_to_search_space()
    - convert_types()
    - get_best_hyperparameters()

Derek van Tilborg | 06-03-2023 | Eindhoven University of Technology

"""

import numpy as np
import os
import torch
from typing import Union
from skopt import gp_minimize
from skopt.space.space import Categorical, Real, Integer
from skopt.utils import use_named_args
from sklearn.metrics import balanced_accuracy_score
from active_learning.hyperparameters import MLP_hypers, GCN_hypers
from active_learning.utils import Evaluate, to_torch_dataloader
from active_learning.acquisition import Acquisition, logits_to_pred
from config import ROOT_DIR


def optimize_hyperparameters(x: np.ndarray, y: np.ndarray, n_calls: int = 50, min_init_points: int = 10,
                             n_folds: int = 5, architecture="gcn", log_file: str = 'hyper_opt.csv') -> dict:
    """ Wrapper function to optimize hyperparameters on a dataset using bootstrapped k-fold cross-validation """

    assert architecture in ['mlp', 'gcn'], f"'architecture' must be 'mlp' or 'gcn'"

    if architecture == 'gcn':
        hypers = GCN_hypers
    if architecture == 'mlp':
        hypers = MLP_hypers

    # Optimize hyperparameters
    opt = BayesianOptimization()
    opt.optimize(x, y, dimensions=hypers, n_calls=n_calls, min_init_points=min_init_points, log_file=log_file,
                 n_folds=n_folds, architecture=architecture)

    best_hypers = get_best_hyperparameters(log_file)

    return best_hypers


class BayesianOptimization:
    def __init__(self):
        """ Init the class with a trainable model. The model class should contain a train() and predict() function
        and be initialized with all of its hyperparameters """
        self.best_score = 1000  # Arbitrary high starting score
        self.history = []
        self.results = None

    def optimize(self, x: np.ndarray, y: np.ndarray, dimensions: dict[str, list[Union[float, str, int]]],
                 n_calls: int = 50, min_init_points: int = 10, log_file: str = 'hypers_search.csv', n_folds: int = 5,
                 architecture: str = 'gcn'):

        # Convert dict of hypers to skopt search_space
        dimensions = {k: [v] if type(v) is not list else v for k, v in dimensions.items()}
        dimensions = dict_to_search_space(dimensions)

        # touch hypers log file
        with open(os.path.join(ROOT_DIR, log_file), 'w') as f:
            f.write(f"score,hypers\n")

        # Objective function for Bayesian optimization
        @use_named_args(dimensions=dimensions)
        def objective(**hyperparameters) -> float:

            # If the same set of hypers gets selected twice (which can happen in the first few runs), skip it
            if hyperparameters in [j for i, j in self.history]:
                score = [i for i, j in self.history if j == hyperparameters][0]
                print(f"skipping - already ran this set of hyperparameters: {hyperparameters}")
            else:
                try:
                    hyperparameters = convert_types(hyperparameters)
                    print(f"Current hyperparameters: {hyperparameters}")

                    score = k_fold_cross_validation(x, y, n_folds=n_folds, architecture=architecture, **hyperparameters)
                    score = 1 - score

                    with open(os.path.join(ROOT_DIR, log_file), 'a') as f:
                        f.write(f"{score},{hyperparameters}\n")

                # If this combination of hyperparameters fails, we use a dummy score that is worse than the best
                except:
                    print(">>  Failed")
                    score = self.best_score + 1

            # append to history and update best score if needed
            self.history.append((score, hyperparameters))
            if score < self.best_score:
                self.best_score = score

            return score

        # Perform Bayesian hyperparameter optimization with n-fold cross-validation
        self.results = gp_minimize(func=objective,
                                   dimensions=dimensions,
                                   acq_func='EI',  # Expected Improvement
                                   n_initial_points=min_init_points,    # Run for this many cycles randomly before BO
                                   n_calls=n_calls,  # Total calls
                                   verbose=True)


def dict_to_search_space(hyperparams: dict[str, list[Union[float, str, int]]]) -> list:
    """ Takes a dict of hyperparameters and converts to skopt search_space"""

    search_space = []

    for k, v in hyperparams.items():
        if type(v[0]) is float and len(v) == 2:
            if k == 'lr' or k == 'learning_rate':
                search_space.append(Real(low=min(v), high=max(v), prior='log-uniform', name=k))
            else:
                search_space.append(Real(low=min(v), high=max(v), name=k))
        elif type(v[0]) is int and len(v) == 2:
            search_space.append(Integer(low=min(v), high=max(v), name=k))
        else:
            search_space.append(Categorical(categories=list(v), name=k))

    return search_space


def convert_types(params: dict) -> dict:
    """ Convert to proper typing. For some reason skopt will mess with float and int typing"""
    new_dict = {}
    for k, v in params.items():
        if isinstance(v, np.generic):
            new_dict[k] = v.item()
        else:
            new_dict[k] = v
    return new_dict


def get_best_hyperparameters(filename: str) -> dict:
    """ Get the best hyperparameters from the log file for an experiment """
    with open(filename) as f:
        best_score = 1000000
        for line in f.readlines()[1:]:
            linesplit = line.split(',')
            if float(linesplit[0]) < best_score:
                hypers_str = ','.join(linesplit[1:]).strip()
                best_score = float(linesplit[0])

    return eval(hypers_str)


def k_fold_cross_validation(x, y, n_folds: int = 5, seed: int = 42, architecture: str = 'mlp', **kwargs) -> float:
    from active_learning.nn import Ensemble

    # Define some variables
    # Set random state and create folds
    rng = np.random.default_rng(seed)
    folds = rng.integers(low=0, high=n_folds, size=len(x))
    scores = []

    for i in range(n_folds):
        # break
        # Subset train/test folds
        if type(x) is np.ndarray:
            x_train, y_train = x[folds != i], y[folds != i]
            x_test, y_test = x[folds == i], y[folds == i]
        else:
            x_train, y_train = [x[j] for j in np.where(folds != i)[0]], y[folds != i]
            x_test, y_test = [x[j] for j in np.where(folds == i)[0]], y[folds == i]

        train_loader = to_torch_dataloader(x_train, y_train, batch_size=256, num_workers=4, shuffle=False)
        test_loader = to_torch_dataloader(x_test, y_test, batch_size=256, num_workers=4, shuffle=False)

        m = Ensemble(architecture=architecture, **kwargs)
        m.train(train_loader, verbose=False)

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
