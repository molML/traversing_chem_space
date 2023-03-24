
import os
from tqdm import tqdm
from active_learning.screening import active_learning
import itertools

LOG_FILE = 'results/uncertainty_based.csv'
PARAMETERS = {'replica': list(range(5)),
              'balanced': [False],
              'max_screen_size': [10000],
              'n_start': [100],
              'batch_size': [32, 64, 128, 256, 512],
              'architecture': ['gcn'],  # , 'bnn'
              'acquisition': ['greedy_exploitative', 'greedy_explorative', 'dynamic_uncertainty']
              }


if __name__ == '__main__':

    experiments = [dict(zip(PARAMETERS.keys(), v)) for v in itertools.product(*PARAMETERS.values())]

    for experiment in tqdm(experiments):
        results = active_learning(n_start=experiment['n_start'],
                                  balanced_start=experiment['balanced'],
                                  acquisition_method=experiment['acquisition'],
                                  max_screen_size=experiment['max_screen_size'],
                                  batch_size=experiment['batch_size'],
                                  architecture=experiment['architecture'],
                                  seed=experiment['replica'])

        # Add the experimental settings to the outfile
        results['acquisition_method'] = experiment['acquisition']
        results['architecture'] = experiment['architecture']
        results['n_start'] = experiment['n_start']
        results['batch_size'] = experiment['batch_size']
        results['replica'] = experiment['replica']
        results['balanced'] = experiment['balanced']

        results.to_csv(LOG_FILE, mode='a', index=False, header=False if os.path.isfile(LOG_FILE) else True)
