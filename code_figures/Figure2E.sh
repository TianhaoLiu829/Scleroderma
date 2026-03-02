#!/bin/bash
#SBATCH --job-name=image_process
#SBATCH --output=image_process_cell_out_smp.log
#SBATCH --error=image_process_cell_error_smp.log
#SBATCH --time=3:00:00             # Adjust based on your job
#SBATCH -c 2
#SBATCH --mem=250G
#SBATCH --cluster=smp

# Load conda and activate environment

module load anaconda3
source activate /ihome/wchen/tianhao/.conda/envs/stardist_env


# Run your Python script
python /ix1/wchen/liutianhao/result/pathology_ST/script/Figure2E.py
