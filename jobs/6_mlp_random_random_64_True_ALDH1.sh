#!/bin/bash
#SBATCH --job-name=6_mlp_random_random_64_True_ALDH1
#SBATCH --output=/home/tilborgd/projects/Active_Learning_Simulation/out/6_mlp_random_random_64_True_ALDH1.out
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=18
#SBATCH --gpus-per-node=1
#SBATCH --time=60:00:00

source $HOME/anaconda3/etc/profile.d/conda.sh
export PYTHONPATH="$PYTHONPATH:$HOME/projects/Active_Learning_Simulation"
$HOME/anaconda3/envs/molml/bin/python -u $HOME/projects/Active_Learning_Simulation/experiments/main.py -o /home/tilborgd/projects/Active_Learning_Simulation/results/6_mlp_random_random_64_True_ALDH1_simulation_results.csv -acq random -bias random -arch mlp -batch_size 64 -retrain True -dataset ALDH1 > $HOME/projects/Active_Learning_Simulation/results/6_mlp_random_random_64_True_ALDH1.log
