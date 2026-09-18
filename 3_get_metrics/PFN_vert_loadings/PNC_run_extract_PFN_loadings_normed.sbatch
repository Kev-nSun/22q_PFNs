#!/bin/bash

#SBATCH --job-name=pfn_extract
#SBATCH --output=PNC_vertex_loading_logs/pfn_extract_%A_%a.out
#SBATCH --error=PNC_vertex_loading_logs/pfn_extract_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-1331%20

set -euo pipefail

mkdir -p PNC_vertex_loading_logs

echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"

python extract_PFN_loadings_normed.py \
    --index "$SLURM_ARRAY_TASK_ID" \
    --csv-path /cbica/projects/bbl_22q/analysis/allometry/inputs/filtered_demo_PNC.csv \
    --norm-dir /cbica/projects/bbl_22q/analysis/allometry/inputs/PNC_PFN_loadings_normed \
    --par-out-dir /cbica/projects/bbl_22q/analysis/pfn/data/PNC_PFNs_nonzero_normed_npz \
    --feat-idx-path /cbica/projects/bbl_22q/analysis/topography/inputs/features_indices.npy \
    --pnc