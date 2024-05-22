import os
import sys
import warnings
import pandas as pd
from active_learning.data_prep import MasterDataset, load_hdf5, get_data, split_data, similarity_vectors
from config import ROOT_DIR
sys.path.append('../active_learning')

warnings.simplefilter(action='ignore', category=FutureWarning)

if __name__ == '__main__':

    # Process the data
    dataset = 'DEMO'

    df = get_data(dataset=dataset)
    df_screen, df_test = split_data(df, screen_size=1000, test_size=200, dataset=dataset)

    MasterDataset(name='screen', df=df_screen, overwrite=True, dataset=dataset)
    MasterDataset(name='test', df=df_test, overwrite=True, dataset=dataset)

    df_screen = pd.read_csv(os.path.join(ROOT_DIR, f'data/{dataset}/original/screen.csv'))
    df_test = pd.read_csv(os.path.join(ROOT_DIR, f'data/{dataset}/original/test.csv'))

    similarity_vectors(df_screen, df_test, dataset=dataset)

