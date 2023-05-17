#!/bin/bash
#SBATCH --job-name=14_gcn_batch_bald_large
#SBATCH --output=/home/tilborgd/projects/Active_Learning_Simulation/out/14_gcn_batch_bald_large.out
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=18
#SBATCH --gpus-per-node=1
#SBATCH --time=120:00:00

source $HOME/anaconda3/etc/profile.d/conda.sh
export PYTHONPATH="$PYTHONPATH:$HOME/projects/Active_Learning_Simulation"
$HOME/anaconda3/envs/molml/bin/python -u $HOME/projects/Active_Learning_Simulation/experiments/main.py -o /home/tilborgd/projects/Active_Learning_Simulation/results -acq batch_bald -bias large -arch gcn > $HOME/projects/Active_Learning_Simulation/results/gcn_batch_bald_large.log
