#!/bin/bash
#SBATCH --job-name=61_gcn_exploration_small_16_True
#SBATCH --output=/home/tilborgd/projects/Active_Learning_Simulation/out/61_gcn_exploration_small_16_True.out
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=18
#SBATCH --gpus-per-node=1
#SBATCH --time=120:00:00

source $HOME/anaconda3/etc/profile.d/conda.sh
export PYTHONPATH="$PYTHONPATH:$HOME/projects/Active_Learning_Simulation"
$HOME/anaconda3/envs/molml/bin/python -u $HOME/projects/Active_Learning_Simulation/experiments/main.py -o /home/tilborgd/projects/Active_Learning_Simulation/results/61_gcn_exploration_small_16_True_simulation_results.csv -acq exploration -bias small -arch gcn -batch_size 16 -retrain True > $HOME/projects/Active_Learning_Simulation/results/61_gcn_exploration_small_16_True.log
