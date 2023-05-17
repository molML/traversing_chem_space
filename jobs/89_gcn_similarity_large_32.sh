#!/bin/bash
#SBATCH --job-name=89_gcn_similarity_large_32
#SBATCH --output=/home/tilborgd/projects/Active_Learning_Simulation/out/89_gcn_similarity_large_32.out
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=18
#SBATCH --gpus-per-node=1
#SBATCH --time=120:00:00

source $HOME/anaconda3/etc/profile.d/conda.sh
export PYTHONPATH="$PYTHONPATH:$HOME/projects/Active_Learning_Simulation"
$HOME/anaconda3/envs/molml/bin/python -u $HOME/projects/Active_Learning_Simulation/experiments/main.py -o /home/tilborgd/projects/Active_Learning_Simulation/results -acq similarity -bias large -arch gcn -batch_size 32 > $HOME/projects/Active_Learning_Simulation/results/gcn_similarity_large_32.log
