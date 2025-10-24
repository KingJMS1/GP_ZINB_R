#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=800M
#SBATCH --time=24:00:00
#SBATCH --output=conv_analysis_model_run_%j_%a.txt
#SBATCH --job-name=conv_analysis_model
#SBATCH --array=1-100

module load anaconda
conda activate renv

srun Rscript model_runner.R ${SLURM_ARRAY_TASK_ID}