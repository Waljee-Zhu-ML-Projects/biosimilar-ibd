#!/bin/bash
#SBATCH --account=awaljee99
#SBATCH --partition=standard
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=8GB
#
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#
#SBATCH -a 1-50

module load Rtidyverse/4.3.1

Rscript 02b_survival_models.R ${SLURM_ARRAY_TASK_ID} ${1} ${2}
