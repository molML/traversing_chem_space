
import numpy as np
from scipy.cluster import hierarchy
import h5py
from typing import Any
import torch
import os
import argparse

ROOT_DIR = os.path.realpath(os.path.dirname(__file__))


def load_hdf5(filename: str) -> Any:
    hf = h5py.File(filename, 'r')
    obj = np.array(hf.get('obj'))
    hf.close()

    return obj


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', help='The path of the output directory', default='results')
    args = parser.parse_args()

    OUT_DIR = args.o

    for dataset in ['PKM2', 'VDR', 'ALDH1']:  # , 'ALDH1'

        D = load_hdf5(f'{OUT_DIR}/tanimoto_distance_vector_{dataset}')

        Z = hierarchy.average(D)

        del D

        torch.save(Z, f'{OUT_DIR}/average_linkage_clustering_{dataset}', pickle_protocol=5)

        del Z

        # / home / tilborgd / projects / random