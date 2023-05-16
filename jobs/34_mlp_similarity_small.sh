#!/bin/bash
#SBATCH --job-name=34_mlp_similarity_small
#SBATCH --output=/home/tilborgd/projects/Active_Learning_Simulation/out/34_mlp_similarity_small.out
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=18
#SBATCH --gpus-per-node=1
#SBATCH --time=100:00:00

source $HOME/anaconda3/etc/profile.d/conda.sh
export PYTHONPATH="$PYTHONPATH:$HOME/projects/Active_Learning_Simulation"
$HOME/anaconda3/envs/molml/bin/python -u $HOME/projects/Active_Learning_Simulation/experiments/snellius_cluster.py -o /home/tilborgd/projects/Active_Learning_Simulation/results -acq similarity -bias small -arch mlp > $HOME/projects/Active_Learning_Simulation/results/mlp_similarity_small.log
