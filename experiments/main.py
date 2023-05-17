
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from tqdm.auto import tqdm
from active_learning.screening import active_learning
import itertools
import argparse


PARAMETERS = {'max_screen_size': [1000],
              'n_start': [64],
              'batch_size': [64, 32, 16],
              'architecture': ['gcn', 'mlp'],
              'seed': list(range(20)),
              'bias': ['random', 'small', 'large'],
              'acquisition': ['random', 'exploration', 'exploitation', 'dynamic', 'batch_bald', 'similarity']
              }


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', help='The path of the output directory', default='results')
    parser.add_argument('-acq', help="Acquisition function ('random', 'exploration', 'exploitation', 'dynamic', "
                                     "'batch_bald', 'similarity')", default='random')
    parser.add_argument('-bias', help='The level of bias ("random", "small", "large")', default='results')
    parser.add_argument('-arch', help='The neural network architecture ("gcn", "mlp")', default='mlp')
    parser.add_argument('-batch_size', help='How many molecules we select each cycle', default=64)
    args = parser.parse_args()

    PARAMETERS['acquisition'] = [args.acq]
    PARAMETERS['bias'] = [args.bias]
    PARAMETERS['architecture'] = [args.arch]
    PARAMETERS['batch_size'] = [args.batch_size]
    LOG_FILE = f'{args.o}/{args.arch}_{args.acq}_{args.bias}_{args.batch_size}_simulation_results.csv'

    # PARAMETERS['acquisition'] = ['random']
    # PARAMETERS['bias'] = ['random']
    # PARAMETERS['architecture'] = ['gcn']
    # LOG_FILE = f"results/{PARAMETERS['architecture'][0]}_{PARAMETERS['acquisition'][0]}_{PARAMETERS['bias'][0]}_simulation_results.csv"

    experiments = [dict(zip(PARAMETERS.keys(), v)) for v in itertools.product(*PARAMETERS.values())]

    for experiment in tqdm(experiments):

        results = active_learning(n_start=experiment['n_start'],
                                  bias=experiment['bias'],
                                  acquisition_method=experiment['acquisition'],
                                  max_screen_size=experiment['max_screen_size'],
                                  batch_size=experiment['batch_size'],
                                  architecture=experiment['architecture'],
                                  seed=experiment['seed'],
                                  optimize_hyperparameters=False)

        # Add the experimental settings to the outfile
        results['acquisition_method'] = experiment['acquisition']
        results['architecture'] = experiment['architecture']
        results['n_start'] = experiment['n_start']
        results['batch_size'] = experiment['batch_size']
        results['seed'] = experiment['seed']
        results['bias'] = experiment['bias']

        results.to_csv(LOG_FILE, mode='a', index=False, header=False if os.path.isfile(LOG_FILE) else True)
