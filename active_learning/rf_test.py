
# Data Processing
import pandas as pd
import numpy as np


from active_learning.data_prep import MasterDataset
from active_learning.data_handler import Handler
from active_learning.utils import Evaluate, to_torch_dataloader
from active_learning.acquisition import Acquisition, logits_to_pred
from active_learning.nn import RfEnsemble
from tqdm.auto import tqdm
from torch.utils.data import DataLoader
from torch.nn import functional as F
from copy import deepcopy

import torch
from torch import Tensor
from torch.utils.data import WeightedRandomSampler
from math import ceil
# Modelling
from sklearn.ensemble import RandomForestClassifier



representation = 'ecfp'
dataset = 'ALDH1'
ds_screen = MasterDataset('screen', representation=representation, dataset=dataset)
ds_test = MasterDataset('test', representation=representation, dataset=dataset)

# logits_N_K_C = [N, num_inference_samples, num_classes]

E = RfEnsemble(10)
E.train(ds_test.x[:1000], ds_test.y[:1000])
logits_N_K_C = E.predict(ds_test.x)

logits_N_K_C.shape


# train_logits_N_K_C = M.predict(train_loader)
# M = Ensemble(seed=seed, ensemble_size=ensemble_size, architecture=architecture, anchored=anchored)
# if cycle == 0 and optimize_hyperparameters:
#     M.optimize_hyperparameters(x_train, y_train)
# M.train(train_loader_balanced, verbose=False)
#
# # Do inference of the train/test/screen data
# print("Train/test/screen inference")
# train_logits_N_K_C = M.predict(train_loader)
# eval_train.eval(train_logits_N_K_C, y_train)











