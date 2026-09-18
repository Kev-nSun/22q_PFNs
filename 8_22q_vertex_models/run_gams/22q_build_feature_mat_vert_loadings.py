#!/usr/bin/env python3
import re
import numpy as np
import pandas as pd
import h5py
from pathlib import Path
"""
This script compiles feature data from npz files outputted fromm "extract_PFN_loadings_normed.py", combining features across subjects, saving it as an h5 file
"""

# Paths
BASE_DIR_PNC = Path("/cbica/projects/bbl_22q/analysis/pfn/data/PNC_PFNs_nonzero_normed_npz")
BASE_DIR_22q = Path("/cbica/projects/bbl_22q/analysis/pfn/data/Final_PFNs_22q_nonzero_normed_npz")
CSV_PATH  = Path("/cbica/projects/bbl_22q/analysis/allometry/inputs/matched_all_demo_QC_median_fd_r1x3.csv")
OUT_H5    = Path("/cbica/projects/bbl_22q/analysis/topography/inputs/22q_PNC_215_matched_compiled_vert_normed_loadings.h5")

# Expected length (from your mask)
TOTAL_LEN = 235355

def normalize_csv_id(x) -> str:
    s = str(x).strip()
    s = s.replace("sub-", "")
    s = re.sub(r"\.0$", "", s)
    if not re.fullmatch(r"\d+", s):
        raise ValueError(f"CSV subject id not purely digits after cleaning: {x} -> {s}")
    return s

def find_subject_dir_by_glob(base_dir: Path, csv_id_digits: str) -> str:
    pattern = f"sub-*{csv_id_digits}"
    candidates = []
    for p in base_dir.glob(pattern):
        if not p.is_dir():
            continue
        if re.fullmatch(r"sub-\d+", p.name) and p.name.split("sub-")[1].endswith(csv_id_digits):
            candidates.append(p.name)

    if len(candidates) == 0:
        raise FileNotFoundError(f"No subject directory found for CSV ID {csv_id_digits} using pattern {pattern}")
    if len(candidates) > 1:
        raise ValueError(f"Multiple subject directories matched CSV ID {csv_id_digits}: {sorted(candidates)}")

    return candidates[0]

def npz_path_for_subject(sub_dirname: str, base_dir: Path) -> Path:
    subj_dir = base_dir / sub_dirname

    # find ses-* directory
    ses_dirs = [p for p in subj_dir.glob("ses-*") if p.is_dir()]
    if len(ses_dirs) == 0:
        raise FileNotFoundError(f"No ses-* directory found under {subj_dir}")
    if len(ses_dirs) > 1:
        raise ValueError(f"Multiple ses-* directories found under {subj_dir}: "
                         f"{[p.name for p in ses_dirs]}")

    ses_dir = ses_dirs[0]
    pfn_dir = ses_dir / "PFN"

    # wildcard the npz filename
    pattern = f"*PFN_loadings_masked_{TOTAL_LEN}.npz"
    matches = list(pfn_dir.glob(pattern))

    if len(matches) == 0:
        raise FileNotFoundError(f"No matching npz found in {pfn_dir} with pattern {pattern}")
    if len(matches) > 1:
        raise ValueError(f"Multiple matching npz files in {pfn_dir}: "
                         f"{[m.name for m in matches]}")

    return matches[0]

def main():
    df = pd.read_csv(CSV_PATH)
    if "subject" not in df.columns:
        raise ValueError("CSV must have a 'subject' column")
    if "stat_22q" not in df.columns:
        raise ValueError("CSV must have a 'stat_22q' column")

    csv_ids = df["subject"].tolist()
    if len(csv_ids) == 0:
        raise ValueError("CSV has no rows")
    
    # Resolve to actual sub-xxxxx directory names in the same order as CSV, BASE_DIR depends on 22q status
    sub_name = []
    for index, id in enumerate(csv_ids):
        digits = normalize_csv_id(id)
        if df["stat_22q"].iloc[index] == 0:
            sub_name.append(find_subject_dir_by_glob(BASE_DIR_PNC, digits))
        else:
            sub_name.append(find_subject_dir_by_glob(BASE_DIR_22q, digits))

    n_subs = len(sub_name)
    print(f"Resolved {n_subs} subject")

    # Make output HDF5
    OUT_H5.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(OUT_H5, "w") as h5:
        # Create dataset: rows=subject, cols=features
        # chunking: chunk along columns for vertex-wise access, rows for writing
        dataset = h5.create_dataset(
            "Y",
            shape=(n_subs, TOTAL_LEN),
            dtype="float32",
            chunks=(1, 8192),      # chunks optimized for writing
            compression="gzip",
            compression_opts=4
        )

        # Save subject ids in the same row order
        h5.create_dataset("subject", data=np.array(sub_name, dtype="S"))

        # Fill dataset
        meta_saved = False

        for index, sub in enumerate(sub_name): #loop through subjects to get npz file (BASE_DIR depends on 22q status)
            if df["stat_22q"].iloc[index] == 0:
                npz_path = npz_path_for_subject(sub, BASE_DIR_PNC)
            else:
                npz_path = npz_path_for_subject(sub, BASE_DIR_22q)
            if not npz_path.exists():
                raise FileNotFoundError(f"Missing npz: {npz_path}")

            data = np.load(npz_path, allow_pickle=False)
            values = data["values"]

            if values.shape[0] != TOTAL_LEN:
                raise ValueError(f"{sub}: values length {values.shape[0]} != {TOTAL_LEN}")

            dataset[index, :] = values.astype(np.float32)

            # Save metadata only once (should be identical across subject)
            if not meta_saved:
                for key in ["features_idx", "net_id", "vertex_id", "vert_num", "total_len"]:
                    if key in data:
                        h5.create_dataset(key, data=data[key])
                meta_saved = True

            if (index + 1) % 50 == 0 or index == n_subs - 1:
                print(f"Wrote {index+1}/{n_subs}")

    print(f"Done. Wrote: {OUT_H5}")

if __name__ == "__main__":
    main()
