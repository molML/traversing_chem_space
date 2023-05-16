#!/bin/bash
#SBATCH --job-name=0_gcn_random_random
#SBATCH --output=/home/tilborgd/projects/Active_Learning_Simulation/out/0_gcn_random_random.out
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=18
#SBATCH --gpus-per-node=1
#SBATCH --time=120:00:00

source $HOME/anaconda3/etc/profile.d/conda.sh
export PYTHONPATH="$PYTHONPATH:$HOME/projects/Active_Learning_Simulation"
$HOME/anaconda3/envs/molml/bin/python -u $HOME/projects/Active_Learning_Simulation/experiments/active_learning.py -o /home/tilborgd/projects/Active_Learning_Simulation/results -acq random -bias random -arch gcn > $HOME/projects/Active_Learning_Simulation/results/gcn_random_random.log