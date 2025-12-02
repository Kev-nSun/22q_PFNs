#!/bin/bash
#SBATCH --job-name=area_calc
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=area_calc_%j.log

export WB_COMMAND=/cbica/software/workbench/bin_linux64/wb_command

python run_area_calcs_PFNs.py