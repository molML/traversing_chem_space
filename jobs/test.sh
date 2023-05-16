#!/bin/bash
#SBATCH --job-name=clustering
#SBATCH --output=/home/tilborgd/projects/Active_Learning_Simulation/out/clustering.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --partition=thin
#SBATCH --time=01:00:00

source $HOME/anaconda3/etc/profile.d/conda.sh
export PYTHONPATH="${PYTHONPATH}:$HOME/projects/Active_Learning_Simulation"
$HOME/anaconda3/envs/molml/bin/python -u $HOME/projects/Active_Learning_Simulation/experiments/snellius_cluster.py -o /home/tilborgd/projects/Active_Learning_Simulation > $HOME/projects/Active_Learning_Simulation/results/clustering.log
