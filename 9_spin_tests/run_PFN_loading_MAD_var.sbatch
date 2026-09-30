#!/bin/bash
#SBATCH --job-name=PFN_MAD
#SBATCH --output=logs/PFN_MAD_%j.out
#SBATCH --error=logs/PFN_MAD_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --tmp=50G

set -euo pipefail

mkdir -p logs

python PFN_loading_MAD_var.py \
    /cbica/projects/bbl_22q/analysis/pfn/data/PFN_loading_MAD_results \
    /cbica/projects/bbl_22q/analysis/pfn/data/PNC_PFNs/Personalized_FN \
    --qc-csv /cbica/projects/bbl_22q/analysis/allometry/inputs/filtered_demo_PNC.csv \
    --chunk-size 4096 \
    --temp-dir "$TMPDIR" \
    --save-npy