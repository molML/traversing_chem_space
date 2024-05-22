import sys
import itertools
import os
import argparse
from tqdm.auto import tqdm
from active_learning.screening import active_learning
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


if __name__ == '__main__':

    PARAMETERS = {'max_screen_size': [128],
                  'n_start': [64],
                  'batch_size': [64],
                  'architecture': ['gcn', 'mlp'],
                  'dataset': ['DEMO'],
                  'seed': [3],
                  'bias': ['random', 'small', 'large'],
                  'acquisition': ['random', 'exploration', 'exploitation', 'dynamic', 'dynamicbald', 'similarity',
                                  'bald']
                  }

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', help='The path of the output directory', default='results')
    parser.add_argument('-acq', help="Acquisition function ('random', 'exploration', 'exploitation', 'dynamic', "
                                     "'similarity')", default='random')
    parser.add_argument('-bias', help='The level of bias ("random", "small", "large")', default='random')
    parser.add_argument('-arch', help='The neural network architecture ("gcn", "mlp")', default='mlp')
    parser.add_argument('-dataset', help='The dataset ("ALDH1", "PKM2", "VDR")', default='DEMO')
    parser.add_argument('-retrain', help='Retrain the model every cycle', default='True')
    parser.add_argument('-batch_size', help='How many molecules we select each cycle', default=64)
    parser.add_argument('-n_start', help='How many molecules we have in our starting set (min=2)', default=64)
    parser.add_argument('-anchored', help='Anchor the weights', default='True')
    args = parser.parse_args()

    PARAMETERS['acquisition'] = [args.acq]
    PARAMETERS['bias'] = ['random']
    PARAMETERS['dataset'] = ['DEMO']
    PARAMETERS['retrain'] = [True]
    PARAMETERS['architecture'] = [args.arch]
    PARAMETERS['batch_size'] = [int(args.batch_size)]
    PARAMETERS['n_start'] = [int(args.n_start)]
    PARAMETERS['anchored'] = [eval(args.anchored)]
    LOG_FILE = args.o

    experiments = [dict(zip(PARAMETERS.keys(), v)) for v in itertools.product(*PARAMETERS.values())]

    for experiment in tqdm(experiments):
        results = active_learning(n_start=experiment['n_start'],
                                  bias=experiment['bias'],
                                  acquisition_method=experiment['acquisition'],
                                  max_screen_size=experiment['max_screen_size'],
                                  batch_size=experiment['batch_size'],
                                  architecture=experiment['architecture'],
                                  seed=experiment['seed'],
                                  retrain=experiment['retrain'],
                                  anchored=experiment['anchored'],
                                  dataset=experiment['dataset'],
                                  optimize_hyperparameters=False)

        # Add the experimental settings to the outfile
        results['acquisition_method'] = experiment['acquisition']
        results['architecture'] = experiment['architecture']
        results['n_start'] = experiment['n_start']
        results['batch_size'] = experiment['batch_size']
        results['seed'] = experiment['seed']
        results['bias'] = experiment['bias']
        results['retrain'] = experiment['retrain']
        results['dataset'] = experiment['dataset']

        results.to_csv(LOG_FILE, mode='a', index=False, header=False if os.path.isfile(LOG_FILE) else True)
