#!/bin/bash
#SBATCH --account=awaljee99
#SBATCH --partition=standard
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=6GB
#
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#
#SBATCH -a 1-50

module load mamba/py3.11
conda activate biosimilar

python 01b_classification_models.py --seed ${SLURM_ARRAY_TASK_ID} --data_file ${1}"_rep"${SLURM_ARRAY_TASK_ID}".csv"