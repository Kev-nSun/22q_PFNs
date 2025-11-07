#!/bin/bash
set -euo pipefail

# Usage:
#   bash cbica/projects/bbl_22q/data/scripts/extract_files_from_datalad.sh <dataset_root> <output_dir> "<file_pattern>"
#
# Example:
#   bash cbica/projects/bbl_22q/data/scripts/extract_files_from_datalad.sh \
#       /cbica/projects/bbl_22q/data/22q_xcpd-0-10-6-BABS_outputs \
#       /cbica/projects/bbl_22q/data/22q_fsLR_91k_dtseries \
#       "xcpd/sub-*/ses-*/func/*_task-idemo_*space-fsLR_den-91k_desc-denoisedSmoothed_bold.dtseries.nii"

ROOT="$1"
DEST="$2"
PATTERN="$3"

mkdir -p "$DEST"

# Sanity: ROOT should be a git/datalad dataset
if ! git -C "$ROOT" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "ERROR: $ROOT is not a git/datalad dataset" >&2
    exit 1
fi

# Prefer annex listing of work-tree files (not .git/annex objects)
if command -v datalad >/dev/null 2>&1; then
    mapfile -t ZIPS < <(git -C "$ROOT" annex find --include '*.zip' --format='${file}\n' | awk 'NF')
else
    mapfile -t ZIPS < <(git -C "$ROOT" ls-files '**/*.zip' | awk 'NF')
fi

if [ "${#ZIPS[@]}" -eq 0 ]; then
    echo "Found 0 zip(s) under $ROOT. Nothing to do."
    exit 0
fi

echo "Found ${#ZIPS[@]} zip(s). Starting extraction…"

for relpath in "${ZIPS[@]}"; do
    # Make absolute path if needed
    case "$relpath" in
        /*) zip_path="$relpath" ;;
        *)  zip_path="$ROOT/$relpath" ;;
    esac

    echo "============================"
    echo "Processing: $zip_path"
    echo "============================"

    # Fetch content (dataset-aware; safe if already present)
    if command -v datalad >/dev/null 2>&1; then
        datalad get -d "$ROOT" "$zip_path" || echo "Note: datalad get not needed for $zip_path"
    fi

    base="$(basename "$zip_path" .zip)"
    out="$DEST/$base"
    mkdir -p "$out"

  # Extract only matching files, preserving archive paths
    if ! 7z x -aos "$zip_path" -o"$out" -ir!"$PATTERN"; then
        echo "Warning: no matches for pattern in $zip_path"
    fi

    # Did anything land in $out?
    if find "$out" -type f | read -r _; then
    # --- Flatten: keep BIDS structure, drop per-zip and any 'xcpd' wrapper ---
        if [ -d "$out/xcpd" ]; then
            echo "Flattening inner xcpd directory from $zip_path ..."
            rsync -a --remove-source-files "$out/xcpd"/ "$DEST"/
        else
            echo "No xcpd/ directory found; copying from $out/"
            rsync -a --remove-source-files "$out"/ "$DEST"/
        fi
    else
        echo "No extracted files in $out (pattern likely matched nothing)."
    fi

    # Clean temp extraction dir
    find "$out" -depth -type d -empty -delete || true
    rmdir "$out" 2>/dev/null || rm -rf "$out"

    # Free space via datalad (dataset-aware). Let DataLad remove local content.
    if command -v datalad >/dev/null 2>&1; then
        datalad drop -d "$ROOT" "$zip_path" || echo "Note: could not drop $zip_path (check remote availability)."
    fi
done

# Optional: prune empties that may remain in DEST after rsync --remove-source-files
find "$DEST" -depth -type d -empty -delete || true

echo "All extractions complete."