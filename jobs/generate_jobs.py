

project_name = "Active_Learning_Simulation"

acquisitions = ['random', 'exploration', 'exploitation', 'bald', 'similarity', 'exploitation_nolearn']  #  'dynamic'
biases = ["random", "small", "large"]
batch_sizes = [64, 32, 16]
architectures = ["gcn", "mlp"]
retrains = ['True', 'False']
datasets = ['ALDH1', 'PKM2', 'VDR']
output = 'results'

# LOG_FILE = f'{args.o}/{args.arch}_{args.acq}_{args.bias}_{args.batch_size}_simulation_results.csv'

# parser = argparse.ArgumentParser()
# parser.add_argument('-o', help='The path of the output directory', default='results')
# parser.add_argument('-acq', help="Acquisition function ('random', 'exploration', 'exploitation', 'dynamic', "
#                                  "'batch_bald', 'similarity')", default='random')
# parser.add_argument('-bias', help='The level of bias ("random", "small", "large")', default='results')
# parser.add_argument('-arch', help='The neural network architecture ("gcn", "mlp")', default='mlp')
# args = parser.parse_args()

experiment = 0
for dataset in datasets:
    for bias in biases:
        for batch_size in batch_sizes:
            if dataset != 'ALDH1' and batch_size != 64:
                continue
            else:

                if batch_size == 64:
                    hours = 16
                elif batch_size == 32:
                    hours = 30
                else:
                    hours = 60

                for arch in architectures:
                    for acq in acquisitions:
                        retrain = 'True'
                        if acq == 'exploitation_nolearn':
                            acq = 'exploitation'
                            retrain = 'False'

                        experiment_name = f"{arch}_{acq}_{bias}_{batch_size}_{retrain}_{dataset}"

                        x = f"""#!/bin/bash
#SBATCH --job-name={experiment}_{experiment_name}
#SBATCH --output=/home/tilborgd/projects/{project_name}/out/{experiment}_{experiment_name}.out
#SBATCH -p gpu
#SBATCH -N 1
#SBATCH --ntasks=18
#SBATCH --gpus-per-node=1
#SBATCH --time={hours}:00:00

source $HOME/anaconda3/etc/profile.d/conda.sh
export PYTHONPATH="${'PYTHONPATH'}:$HOME/projects/{project_name}"
$HOME/anaconda3/envs/molml/bin/python -u $HOME/projects/{project_name}/experiments/main.py -o /home/tilborgd/projects/{project_name}/{output}/{experiment}_{experiment_name}_simulation_results.csv -acq {acq} -bias {bias} -arch {arch} -batch_size {batch_size} -retrain {retrain} -dataset {dataset} > $HOME/projects/{project_name}/{output}/{experiment}_{experiment_name}.log
"""


                        filename = f"jobs/{experiment}_{experiment_name}.sh"
                        print(filename)

                        # write x to file
                        with open(filename, 'w') as f:
                            f.write(x)

                        experiment += 1

