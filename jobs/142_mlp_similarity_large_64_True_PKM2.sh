#!/bin/bash
#SBATCH --job-name=142_mlp_similarity_large_64_True_PKM2
#SBATCH --output=/home/tilborgd/projects/Active_Learning_Simulation/out/142_mlp_similarity_large_64_True_PKM2.out
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=18
#SBATCH --gpus-per-node=1
#SBATCH --time=60:00:00

source $HOME/anaconda3/etc/profile.d/conda.sh
export PYTHONPATH="$PYTHONPATH:$HOME/projects/Active_Learning_Simulation"
$HOME/anaconda3/envs/molml/bin/python -u $HOME/projects/Active_Learning_Simulation/experiments/main.py -o /home/tilborgd/projects/Active_Learning_Simulation/results/142_mlp_similarity_large_64_True_PKM2_simulation_results.csv -acq similarity -bias large -arch mlp -batch_size 64 -retrain True -dataset PKM2 > $HOME/projects/Active_Learning_Simulation/results/142_mlp_similarity_large_64_True_PKM2.log
