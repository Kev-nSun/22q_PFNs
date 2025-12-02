#!/usr/bin/env bash
set -e

# Subject+session loop for T1 QC PNGs on CUBIC (using raw BIDS T1)
BIDS="/cbica/projects/bbl_22q/data/raw_datalad"
SIF="/cbica/projects/bbl_22q/babs/qsiprep-1.0.1.sif"
OUT="/cbica/projects/bbl_22q/data/T1_QC/slices"
AFF="/cbica/projects/bbl_22q/data/T1_QC/input_files"
TPL="/cbica/projects/bbl_22q/templateflow_home/tpl-MNI152NLin6Asym/tpl-MNI152NLin6Asym_res-01_T1w.nii.gz"

mkdir -p "$OUT" "$AFF"

for h5 in $(ls /cbica/projects/bbl_22q/data/derivatives_2025/22q_fmriprep_anat_extracted/fmriprep_anat/sub-*/ses-*/anat/*from-T1w_to-MNI152NLin6Asym*_xfm.h5); do
  anatdir=$(dirname "$h5")
  sessdir=$(dirname "$anatdir")
  subdir=$(dirname "$sessdir")
  ses=$(basename "$sessdir")
  sub=$(basename "$subdir")
  t1="$BIDS/$sub/$ses/anat/${sub}_${ses}_T1w.nii.gz"
  t1=$(ls $t1 | head -n1)
  # Ensure the raw T1 is materialized if BIDS is a DataLad dataset
  if command -v datalad >/dev/null 2>&1; then datalad -C "$BIDS" get "$t1"; fi

  # Ensure the raw T1 is present from DataLad
  datalad -C "$BIDS" get "$t1"

  # Disassemble in the anat dir (either --pwd or bash -lc)
  apptainer exec -e -B /cbica/projects/bbl_22q/ --pwd "$AFF" "$SIF" \
  CompositeTransformUtil --disassemble "$h5" "disassembled_${sub}_${ses}_"
  # Pick the affine from that anat dir (index like 00_ is now handled by the glob)
  affine=$(ls "$AFF"/*disassembled_${sub}_${ses}_*Affine*.mat 2>/dev/null | head -n1)

  outnii="${OUT}/${sub}_${ses}_space-MNI152NLin6Asym_T1w.nii.gz"
  apptainer exec -e -B /cbica/projects/bbl_22q/ "$SIF" \
  antsApplyTransforms -d 3 -i "$t1" -t "$affine" -r "$TPL" -o "$outnii" --interpolation lanczoswindowedsinc

  slicer "$outnii" \
    -x 0.4 "${OUT}/${sub}_${ses}_S1.png" \
    -x 0.6 "${OUT}/${sub}_${ses}_S3.png" \
    -z 0.5 "${OUT}/${sub}_${ses}_A2.png" \
    -z 0.6 "${OUT}/${sub}_${ses}_A3.png"
done
