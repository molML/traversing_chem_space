

![python version](https://img.shields.io/badge/python-v.3.9-blue)
![license](https://img.shields.io/badge/license-MIT-orange)


<h1 id="benchmark-study">Traversing Chemical Space with Active Deep Learning</h1>
Derek van Tilborg and Francesca Grisoni*

*f.grisoni@tue.nl

<h2 id="benchmark-study">Abstract</h2>

Deep learning is accelerating drug discovery. However, current approaches are often affected by limitations in the available data, e.g., in terms of size or molecular diversity. Active deep learning has
an untapped potential for low-data drug discovery, as it allows to improve a model iteratively during the screening process by acquiring new data, and to adjust its course along the way. However, several known unknowns exist when it comes to active learning: (a) what the best computational strategies are for chemical space exploration, (b) how active learning holds up to traditional, non-iterative, approaches, and (c) how it should be used in the low-data scenarios typical of drug discovery. These open questions currently limit the wider adoption of active learning in drug discovery. To provide answers, this study simulates a real-world low-data drug discovery scenario, and systematically analyses six active learning strategies combined with two deep learning architectures, on three large- scale molecular libraries. Not only do we show that active learning can achieve up to a six-fold improvement in hit discovery compared to traditional methods, but we also identify the most important determinants of its success in low-data regimes. This study lays the first-in-time foundations for the prospective use of active deep learning for low-data drug discovery and is expected to accelerate its adoption.



![Figure 1](figures/fig1.png)

 
## Description
This repository contains all code to replicate our study. Feel free to try out anything you find here.

## Modules
- `data_prep.py`: Processes active and inactive compound data.
- `clustering.py`: Clusters compounds for sampling diversity.
- `nn.py`: Contains neural network models (MLP, GCN, etc.).
- `screening.py`: Core script for active learning cycles.
- `utils.py`: Utility functions for data handling and evaluation.
- `main.py`: Entry point for running experiments with customizable parameters.
 
## Usage
Run `python main.py` with desired command-line arguments to start the active learning process. Ensure necessary data files (active and inactive `.smi`) are present.
 
## Requirements
Install dependencies from the provided env.yaml file.

```conda env create -f env.yaml```

This codebase uses Python 3.9 and depends on:
- [PyTorch](https://pytorch.org/) (2.0.1)
- [PyTorch Geometric](https://pytorch-geometric.readthedocs.io/en/latest/) (2.3.1)
- [RDKit](https://www.rdkit.org/) (2023.3.2)
- [Scikit-learn](https://scikit-learn.org/) (1.3.0)


<!-- How to cite-->
<h2 id="How-to-cite">How to cite</h2>
You can currently cite our preprint

Traversing Chemical Space with Active Deep Learning. Derek van Tilborg and Francesca Grisoni.
ChemRxiv, 2023.
DOI: https://doi.org/10.26434/chemrxiv-2023-wgl32-v3


<!-- License-->
<h2 id="License">License</h2>

All code is under MIT license.
