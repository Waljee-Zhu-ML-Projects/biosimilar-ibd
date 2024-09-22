#!/bin/bash
#SBATCH --account=awaljee99
#SBATCH --partition=standard
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=16GB
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load Rtidyverse/4.3.1

Rscript ${1}.R ${2} ${3} ${4} ${5} ${6}