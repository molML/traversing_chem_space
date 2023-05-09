""" Range of hyperparameters that will be optimized for each model"""

MLP_hypers = {'lr': [1e-3, 1e-4, 1e-5],                 # categorical
              'n_hidden': [1024],                       # categorical
              'epochs': [50, 100, 250, 500],
              'ensemble_size': [10]}                    # categorical


GCN_hypers = {'lr': [1e-3, 1e-4, 1e-5],                 # categorical
              'n_hidden': [1024],                       # categorical
              'epochs': [50, 100, 250],                 # categorical
              'num_conv_layers': [3, 5],
              'ensemble_size': [10]}                    # integer range from 2 - 5

