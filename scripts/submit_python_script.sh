#!/bin/bash
#SBATCH --account=awaljee99
#SBATCH --partition=standard
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=16GB
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load mamba/py3.11
conda activate biosimilar

python ${1}.py
