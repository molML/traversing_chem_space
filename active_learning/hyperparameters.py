""" Range of hyperparameters that are optimized for each model"""

BNN_hypers = {'lr': [1e-3, 1e-4, 1e-5],                 # categorical
              'hidden_size': [128, 256, 512, 1024],     # categorical
              'epochs': [100, 1000, 5000, 10000]}       # categorical


GCN_hypers = {'lr': [1e-3, 1e-4, 1e-5],                 # categorical
              'hidden_channels': [128, 256, 512, 1024],     # categorical
              'epochs': [50, 100, 200, 500],            # categorical
              'num_conv_layers': [1, 5]}                # integer range from 2 - 5


BGCN_hypers = {}
