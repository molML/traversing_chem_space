""" Range of hyperparameters that will be optimized for each model"""

MLP_hypers = {'lr': [3e-3, 3e-4, 3e-5],                 # categorical
              'n_hidden': [512, 1024],                  # categorical
              'epochs': [100],
              'anchored': [True],
              'ensemble_size': [5]}                     # categorical


GCN_hypers = {'lr': [3e-3, 3e-4, 3e-5],                 # categorical
              'n_hidden': [512, 1024],                  # categorical
              'epochs': [100],                          # categorical
              'anchored': [True],
              'num_conv_layers': [3],
              'ensemble_size': [5]}                     # integer range from 2 - 5
