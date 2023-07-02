#!/bin/bash
#SBATCH --job-name=15_gcn_bald_random_32_True
#SBATCH --output=/home/tilborgd/projects/Active_Learning_Simulation/out/15_gcn_bald_random_32_True.out
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=18
#SBATCH --gpus-per-node=1
#SBATCH --time=120:00:00

source $HOME/anaconda3/etc/profile.d/conda.sh
export PYTHONPATH="$PYTHONPATH:$HOME/projects/Active_Learning_Simulation"
$HOME/anaconda3/envs/molml/bin/python -u $HOME/projects/Active_Learning_Simulation/experiments/main.py -o /home/tilborgd/projects/Active_Learning_Simulation/results/15_gcn_bald_random_32_True_simulation_results.csv -acq bald -bias random -arch gcn -batch_size 32 -retrain True > $HOME/projects/Active_Learning_Simulation/results/15_gcn_bald_random_32_True.log
