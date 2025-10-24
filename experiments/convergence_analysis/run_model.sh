#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --nodes=25
#SBATCH --ntasks=100
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=conv_analysis_model_run_%j.txt
#SBATCH --job-name=conv_analysis_model

module load anaconda
conda activate renv

srun --multi-prog run_model.conf