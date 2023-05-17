#!/bin/bash
#SBATCH --job-name=60_gcn_random_small_16
#SBATCH --output=/home/tilborgd/projects/Active_Learning_Simulation/out/60_gcn_random_small_16.out
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=18
#SBATCH --gpus-per-node=1
#SBATCH --time=120:00:00

source $HOME/anaconda3/etc/profile.d/conda.sh
export PYTHONPATH="$PYTHONPATH:$HOME/projects/Active_Learning_Simulation"
$HOME/anaconda3/envs/molml/bin/python -u $HOME/projects/Active_Learning_Simulation/experiments/main.py -o /home/tilborgd/projects/Active_Learning_Simulation/results -acq random -bias small -arch gcn -batch_size 16 > $HOME/projects/Active_Learning_Simulation/results/gcn_random_small_16.log
