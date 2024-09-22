#!/bin/bash
#SBATCH --account=awaljee99
#SBATCH --partition=standard
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=1GB
#
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#
#SBATCH -a 1-50

module load Rtidyverse/4.3.1

Rscript 01a_classification_data.R ${SLURM_ARRAY_TASK_ID} ${1} ${2} ${3}
