#!/bin/bash
#PBS -N scRNA_noise
#PBS -q x86_64
#PBS -l nodes=1:ppn=16,mem=128gb,walltime=06:00:00
#PBS -o /home/mllorens/mbeckel/transcriptional_noise/logs/noise_analysis.log
#PBS -m abe
#PBS -M ms.beckel@gmail.com

# Activate Conda
source /ngs/mllorens/mbeckel/miniconda3/bin/activate
conda activate r_env

# Run script
Rscript /home/mllorens/mbeckel/projects/transcriptional_noise/scripts/buckley2022_tn_analysis.R