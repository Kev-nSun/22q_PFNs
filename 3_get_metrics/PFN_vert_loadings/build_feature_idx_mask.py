#!/usr/bin/env python3

import argparse
import re
import numpy as np
import pandas as pd
import scipy.io
from pathlib import Path


"""
Can be called using these arguments:

python build_feature_idx_mask.py \
    --csv-path /cbica/projects/bbl_22q/analysis/allometry/inputs/filtered_demo_PNC.csv \
    --fn-dir /cbica/projects/bbl_22q/analysis/pfn/data/PNC_PFNs/Personalized_FN \
    --out-path /cbica/projects/bbl_22q/analysis/allometry/inputs/PNC_features_indices.npy

Builds a PFN feature mask using PNC subjects.

For every PNC subject in filtered_demo_PNC.csv:
    1. Load the subject's raw FN.mat
    2. Ensure FN has shape (59412 vertices, 17 PFNs)
    3. Flatten in network-major order:
           PFN01 vertex 1 ... vertex 59412
           PFN02 vertex 1 ... vertex 59412
           ...
           PFN17 vertex 1 ... vertex 59412
    4. Retain any feature that is nonzero in at least one PNC subject.

Output:
    features_indices.npy

The saved indices are 0-based indices into the full flattened
17 * 59412 = 1,010,004 dimensional PFN feature space.

This mask can then be supplied to the PFN extraction script.

IMPORTANT:
    The mask is intentionally based on raw FN.mat files.
    Vertex-wise normalization does not affect whether an individual
    loading is zero, so normalization is not required for this step.
"""


VERT_NUM = 59412
N_NETWORKS = 17
N_FEATURES = VERT_NUM * N_NETWORKS


def normalize_csv_id(x) -> str:
    """
    Normalize participant IDs to digits only.

    Examples:
        sub-1234 -> 1234
        1234     -> 1234
        1234.0   -> 1234
    """
    s = str(x).strip()
    s = s.replace("sub-", "")
    s = re.sub(r"\.0$", "", s)

    if not re.fullmatch(r"\d+", s):
        raise ValueError(
            f"CSV subject ID not purely digits after cleaning: "
            f"{x} -> {s}"
        )

    return s


def find_subject_dir_by_glob(
    fn_dir: Path,
    csv_id_digits: str,
) -> Path:
    """
    Find a unique PNC subject directory under fn_dir.

    Allows the CSV ID to match the end of the directory ID.

    Example:
        CSV ID: 992
        directory: sub-0992
    """

    pattern = f"sub-*{csv_id_digits}"
    candidates = []

    for p in fn_dir.glob(pattern):

        if not p.is_dir():
            continue

        if not re.fullmatch(r"sub-\d+", p.name):
            continue

        subject_digits = p.name.replace("sub-", "")

        if subject_digits.endswith(csv_id_digits):
            candidates.append(p)

    if len(candidates) == 0:
        raise FileNotFoundError(
            f"No subject directory found for CSV ID "
            f"{csv_id_digits} using pattern {pattern}"
        )

    if len(candidates) > 1:
        raise ValueError(
            f"Multiple subject directories matched CSV ID "
            f"{csv_id_digits}: "
            f"{sorted(p.name for p in candidates)}"
        )

    return candidates[0]


def choose_fn_array_from_mat(matdict: dict) -> np.ndarray:
    """
    Select the PFN matrix from FN.mat.

    Preference:
        1. MATLAB variable named 'FN'
        2. Otherwise, a unique numeric array with shape:
               (59412, 17)
           or
               (17, 59412)
    """

    if "FN" in matdict:
        arr = matdict["FN"]

        if isinstance(arr, np.ndarray):
            return arr

    candidates = []

    for key, value in matdict.items():

        if key.startswith("__"):
            continue

        if not isinstance(value, np.ndarray):
            continue

        if value.ndim != 2:
            continue

        if not np.issubdtype(value.dtype, np.number):
            continue

        if (
            value.shape == (VERT_NUM, N_NETWORKS)
            or value.shape == (N_NETWORKS, VERT_NUM)
        ):
            candidates.append((key, value))

    if len(candidates) == 1:
        return candidates[0][1]

    if len(candidates) > 1:
        keys = [key for key, _ in candidates]

        raise ValueError(
            f"Multiple plausible FN arrays found: {keys}"
        )

    raise ValueError(
        "Could not find an FN array with shape "
        f"({VERT_NUM}, {N_NETWORKS}) or "
        f"({N_NETWORKS}, {VERT_NUM})."
    )


def load_fn_mat(fn_path: Path) -> np.ndarray:
    """
    Load FN.mat and return FN with shape:

        vertices x networks
        (59412, 17)
    """

    if not fn_path.exists():
        raise FileNotFoundError(
            f"FN.mat not found: {fn_path}"
        )

    try:
        mat = scipy.io.loadmat(str(fn_path))

    except NotImplementedError as e:
        raise NotImplementedError(
            f"{fn_path} appears to be MATLAB v7.3/HDF5. "
            "scipy.io.loadmat cannot read it."
        ) from e

    arr = np.asarray(
        choose_fn_array_from_mat(mat)
    )

    if arr.shape == (VERT_NUM, N_NETWORKS):
        fn = arr

    elif arr.shape == (N_NETWORKS, VERT_NUM):
        fn = arr.T

    else:
        raise ValueError(
            f"Unexpected FN shape {arr.shape} in {fn_path}. "
            f"Expected ({VERT_NUM}, {N_NETWORKS}) "
            f"or ({N_NETWORKS}, {VERT_NUM})."
        )

    return fn.astype(np.float32, copy=False)


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Build features_indices.npy from raw PNC FN.mat files."
        )
    )

    parser.add_argument(
        "--csv-path",
        type=str,
        default=(
            "/cbica/projects/bbl_22q/analysis/allometry/"
            "inputs/filtered_demo_PNC.csv"
        ),
        help=(
            "PNC phenotype CSV containing participant_id."
        ),
    )

    parser.add_argument(
        "--fn-dir",
        type=str,
        default=(
            "/cbica/projects/bbl_22q/analysis/pfn/data/"
            "PNC_PFNs/Personalized_FN"
        ),
        help=(
            "PNC PFN directory containing "
            "sub-XXXX/FN.mat."
        ),
    )

    parser.add_argument(
        "--out-path",
        type=str,
        default=(
            "/cbica/projects/bbl_22q/analysis/topography/"
            "inputs/features_indices.npy"
        ),
        help="Output path for features_indices.npy.",
    )

    args = parser.parse_args()

    csv_path = Path(args.csv_path)
    fn_dir = Path(args.fn_dir)
    out_path = Path(args.out_path)

    # ---------------------------------------------------------
    # Validate input paths
    # ---------------------------------------------------------

    if not csv_path.exists():
        raise FileNotFoundError(
            f"PNC CSV not found: {csv_path}"
        )

    if not fn_dir.exists():
        raise FileNotFoundError(
            f"PNC PFN directory not found: {fn_dir}"
        )

    # ---------------------------------------------------------
    # Load PNC subjects
    # ---------------------------------------------------------

    df = pd.read_csv(csv_path)

    if "participant_id" not in df.columns:
        raise ValueError(
            "PNC CSV must contain a 'participant_id' column."
        )

    csv_subjects = df["participant_id"].tolist()

    if len(csv_subjects) == 0:
        raise ValueError(
            f"No PNC participants found in {csv_path}"
        )

    print(f"PNC subjects in CSV: {len(csv_subjects):,}")
    print(f"Vertices per PFN:     {VERT_NUM:,}")
    print(f"Number of PFNs:       {N_NETWORKS}")
    print(f"Full feature count:   {N_FEATURES:,}")

    # ---------------------------------------------------------
    # Initialize cohort-wide mask
    #
    # False:
    #     feature has been zero for every subject seen so far
    #
    # True:
    #     feature is nonzero in >= 1 PNC subject
    # ---------------------------------------------------------

    nonzero_mask = np.zeros(
        N_FEATURES,
        dtype=bool,
    )

    valid_subjects = []
    missing_subjects = []

    # ---------------------------------------------------------
    # Scan raw FN.mat files
    # ---------------------------------------------------------

    for index, csv_subject in enumerate(csv_subjects):

        subject_digits = normalize_csv_id(
            csv_subject
        )

        try:

            subject_dir = find_subject_dir_by_glob(
                fn_dir,
                subject_digits,
            )

        except FileNotFoundError:

            print(
                f"WARNING: subject directory not found for "
                f"{csv_subject}; skipping"
            )

            missing_subjects.append(
                str(csv_subject)
            )

            continue

        fn_path = subject_dir / "FN.mat"

        if not fn_path.exists():

            print(
                f"WARNING: FN.mat missing for "
                f"{subject_dir.name}; skipping"
            )

            missing_subjects.append(
                str(csv_subject)
            )

            continue

        # Load raw PFN loadings:
        #
        # shape = 59412 vertices x 17 PFNs
        fn = load_fn_mat(fn_path)

        # -----------------------------------------------------
        # Flatten into network-major feature ordering.
        #
        # fn.T:
        #       17 x 59412
        #
        # ravel():
        #
        # index 0                  -> PFN01 vertex 0
        # index 59411             -> PFN01 vertex 59411
        # index 59412             -> PFN02 vertex 0
        # ...
        # index 1,010,003         -> PFN17 vertex 59411
        #
        # This matches the extraction script's:
        #
        # network = feature_idx // 59412
        # vertex  = feature_idx % 59412
        # -----------------------------------------------------

        flat = fn.T.ravel()

        if flat.size != N_FEATURES:
            raise ValueError(
                f"{subject_dir.name}: flattened FN contains "
                f"{flat.size:,} features; expected "
                f"{N_FEATURES:,}"
            )

        # Retain the feature if it is nonzero for this
        # participant OR any previously processed participant.
        #
        # After all subjects this is equivalent to:
        #
        # ~np.all(features == 0, axis=0)

        nonzero_mask |= (flat != 0)

        valid_subjects.append(
            subject_dir.name
        )

        if (
            (index + 1) % 50 == 0
            or index == len(csv_subjects) - 1
        ):

            print(
                f"Scanned {index + 1:,}/"
                f"{len(csv_subjects):,} CSV subjects | "
                f"valid={len(valid_subjects):,} | "
                f"features retained="
                f"{np.count_nonzero(nonzero_mask):,}"
            )

    # ---------------------------------------------------------
    # Validate cohort
    # ---------------------------------------------------------

    if len(valid_subjects) == 0:
        raise RuntimeError(
            "No valid PNC FN.mat files were found. "
            "Mask was not created."
        )

    # ---------------------------------------------------------
    # Convert Boolean mask into original feature indices
    # ---------------------------------------------------------

    kept_feature_indices = np.flatnonzero(
        nonzero_mask
    ).astype(np.int32)

    # ---------------------------------------------------------
    # Save features_indices.npy
    # ---------------------------------------------------------

    out_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    np.save(
        out_path,
        kept_feature_indices,
        allow_pickle=False,
    )

    # ---------------------------------------------------------
    # Summary
    # ---------------------------------------------------------

    n_kept = kept_feature_indices.size
    n_removed = N_FEATURES - n_kept

    print()
    print("========================================")
    print("PNC PFN feature mask complete")
    print("========================================")
    print(
        f"PNC subjects in CSV:   "
        f"{len(csv_subjects):,}"
    )
    print(
        f"Valid FN.mat files:    "
        f"{len(valid_subjects):,}"
    )
    print(
        f"Missing/skipped:       "
        f"{len(missing_subjects):,}"
    )
    print(
        f"Full feature count:    "
        f"{N_FEATURES:,}"
    )
    print(
        f"Retained features:     "
        f"{n_kept:,}"
    )
    print(
        f"All-zero features:     "
        f"{n_removed:,}"
    )
    print(
        f"Percent retained:      "
        f"{100 * n_kept / N_FEATURES:.2f}%"
    )
    print(
        f"Saved:                 "
        f"{out_path}"
    )


if __name__ == "__main__":
    main()