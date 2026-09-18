#!/bin/bash
#SBATCH --job-name=temp_area_calc_22q
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=temp_area_calc_22q_%j.log

export WB_COMMAND=/cbica/software/workbench/bin_linux64/wb_command

python 22q_run_template_area_calcs_PFNs.py