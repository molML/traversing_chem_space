
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from tqdm.auto import tqdm
from active_learning.screening import active_learning
import itertools


PARAMETERS = {'max_screen_size': [1000],
              'n_start': [64],
              'batch_size': [64],
              'architecture': ['rf'],
              'dataset': ['ALDH1', 'PKM2', 'VDR', 'IDH1', 'ADRB2', 'OPRK1', 'GBA', 'KAT2A', 'FEN1'],
              'bias': ['random'],
              'retrain': [True],
              'scrambledx': [False],
              'scrambledx_seed': [1, 2, 3],
              'acquisition': ['random', 'exploration', 'exploitation', 'similarity', 'bald'],
              }


if __name__ == '__main__':

    # Setup the grid of experiments
    experiments_ = [dict(zip(PARAMETERS.keys(), v)) for v in itertools.product(*PARAMETERS.values())]
    experiments = []
    for exp in experiments_:
        if exp['retrain'] == False:
            if exp['acquisition'] == 'exploitation':
                experiments.append(exp)
        else:
            experiments.append(exp)

    # for i, experiment in tqdm(zip(list(range(len(experiments))), experiments)):
        # print(i, experiment)
    for i, experiment in tqdm(zip(list(range(len(experiments)))[::-1], experiments[::-1])):
        print(i, experiment)
        LOG_FILE = f"results/{i}_{experiment['architecture']}_{experiment['acquisition']}_{experiment['bias']}_{experiment['n_start']}_{experiment['retrain']}_{experiment['dataset']}_shuffle_{experiment['scrambledx_seed']}_simulation_results.csv"
        # LOG_FILE = f"results/{i}_{experiment['architecture']}_{experiment['acquisition']}_{experiment['bias']}_{experiment['n_start']}_{experiment['retrain']}_{experiment['dataset']}_logfixed_simulation_results.csv"

        if os.path.exists(LOG_FILE):
                print(f"Experiment {i} already exists")
        else:
            for seed in range(10):

                results = active_learning(n_start=experiment['n_start'],
                                          bias=experiment['bias'],
                                          acquisition_method=experiment['acquisition'],
                                          max_screen_size=experiment['max_screen_size'],
                                          batch_size=experiment['batch_size'],
                                          architecture=experiment['architecture'],
                                          seed=seed,
                                          retrain=experiment['retrain'],
                                          dataset=experiment['dataset'],
                                          scrambledx=experiment['scrambledx'],
                                          scrambledx_seed=experiment['scrambledx_seed'],
                                          optimize_hyperparameters=False)

                # Add the experimental settings to the outfile
                results['acquisition_method'] = experiment['acquisition']
                results['architecture'] = experiment['architecture']
                results['n_start'] = experiment['n_start']
                results['batch_size'] = experiment['batch_size']
                results['seed'] = seed
                results['bias'] = experiment['bias']
                results['retrain'] = experiment['retrain']
                results['scrambledx'] = experiment['scrambledx']
                results['scrambledx_seed'] = experiment['scrambledx_seed']
                results['dataset'] = experiment['dataset']

                results.to_csv(LOG_FILE, mode='a', index=False, header=False if os.path.isfile(LOG_FILE) else True)
