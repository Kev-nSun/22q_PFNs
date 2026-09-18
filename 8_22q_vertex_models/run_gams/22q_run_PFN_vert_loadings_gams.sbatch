#!/bin/bash
#SBATCH --job-name=pfn_vert_gams
#SBATCH --output=vert_loading_gam_logTC_logs/pfn_gams_%A_%a.out
#SBATCH --error=vert_loading_gam_logTC_logs/pfn_gams_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-235%50

set -euo pipefail

module load neuroR/0.2.0

H5="/cbica/projects/bbl_22q/analysis/topography/inputs/22q_PNC_215_matched_compiled_vert_normed_loadings.h5"
CSV="/cbica/projects/bbl_22q/analysis/allometry/inputs/22q_PNC_matched_all.csv"
OUTDIR="/cbica/projects/bbl_22q/analysis/topography/results/091426_normed_22q_vert_gam_results_nonstd_sex_avg_logTC_subclass"

mkdir -p vert_loading_gam_logTC_logs "${OUTDIR}"

DONEFILE="${OUTDIR}/gam_chunk_$(printf "%04d" ${SLURM_ARRAY_TASK_ID})_DONE.txt"
if [ -f "${DONEFILE}" ]; then
  echo "Chunk ${SLURM_ARRAY_TASK_ID}: DONE exists, skipping."
  exit 0
fi

echo "Job: ${SLURM_JOB_ID}  Task: ${SLURM_ARRAY_TASK_ID}"
echo "Host: $(hostname)"
echo "H5: ${H5}"
echo "CSV: ${CSV}"
echo "OUTDIR: ${OUTDIR}"

Rscript 22q_PFN_vert_loadings_gams_chunk.R \
  --h5 "${H5}" \
  --csv "${CSV}" \
  --outdir "${OUTDIR}" \
  --chunk-id "${SLURM_ARRAY_TASK_ID}" \
  --chunk-size 1000

touch "${DONEFILE}"