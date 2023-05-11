#!/bin/bash
#SBATCH --job-name=clustering
#SBATCH --output=/home/tilborgd/projects/random/out/clustering.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --partition=thin
#SBATCH --time=01:00:00

source $HOME/anaconda3/etc/profile.d/conda.sh
export PYTHONPATH="${PYTHONPATH}:$HOME/projects/random"
$HOME/anaconda3/envs/molml/bin/python -u $HOME/projects/random/snellius_cluster.py -o /home/tilborgd/projects/random > $HOME/projects/random/clustering.log