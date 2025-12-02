#!/bin/bash
#SBATCH --job-name=PFN_norm
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=/cbica/projects/bbl_22q/analysis/allometry/scripts/PNC_norm_loadings_logs/PFN_norm_%A_%a.log
#SBATCH --array=1-1482

LOG_DIR=/cbica/projects/bbl_22q/analysis/allometry/scripts/PNC_norm_loadings_logs
mkdir -p "$LOG_DIR"

set -euo pipefail
module load neuroR/0.2.0
which Rscript

echo "JOB $SLURM_JOB_ID TASK $SLURM_ARRAY_TASK_ID"
echo "Rscript Norm_PFNs_sum_to_1.R --array_id $SLURM_ARRAY_TASK_ID"

Rscript /cbica/projects/bbl_22q/analysis/allometry/scripts/Norm_PFNs_sum_to_1.R \
    --subject_list /cbica/projects/bbl_22q/analysis/allometry/inputs/PNC_subject_list.txt \
    --in_dir /cbica/projects/bbl_22q/analysis/pfn/data/PNC_PFNs/Personalized_FN \
    --out_dir /cbica/projects/bbl_22q/analysis/allometry/inputs/PNC_PFN_loadings_normed \
    --array_id "$SLURM_ARRAY_TASK_ID"