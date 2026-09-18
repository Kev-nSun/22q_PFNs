#!/usr/bin/env python3
import argparse
import re
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Tuple

"""
Extract masked, vertex-wise normalized PFN loadings from normalized CSV files.

The normalized CSVs are produced by Norm_PFNs_sum_to_1_and_write_ciftis.R and
are expected at:

    <norm-dir>/<sub_name>/<sub_name>_PFN_loadings_normed.csv

Each CSV must contain 59,412 rows (vertices) and 17 columns (PFNs).

A precomputed PNC-derived features_indices.npy defines which entries of the
full 17 x 59,412 PFN feature space are retained. The same mask can be applied
to both PNC and 22q subjects.

Feature ordering is network-major:
    PFN01 vertices 0..59411,
    PFN02 vertices 0..59411,
    ...
    PFN17 vertices 0..59411.

Outputs remain one compressed .npz per subject/session.
"""

VERT_NUM = 59412
N_NETWORKS = 17


def normalize_csv_id(x) -> str:
    s = str(x).strip()
    s = s.replace("sub-", "").replace("ses-", "")
    s = re.sub(r"\.0$", "", s)
    if not re.fullmatch(r"\d+", s):
        raise ValueError(
            f"CSV subject/session id not purely digits after cleaning: {x} -> {s}"
        )
    return s


def find_subject_dir_by_glob(
    base_area_dir: Path,
    csv_id_digits: str,
    csv_ses_digits: str,
) -> Tuple[str, str]:
    """
    Resolve a 22q CSV subject/session to the directory names used downstream.

    Finds:
        sub-<digits ending with csv_id_digits>/
            ses-<digits ending with csv_ses_digits>

    Returns:
        (sub_name, ses_name)
    """
    subj_pattern = f"sub-*{csv_id_digits}"
    subject_candidates = []

    for p in base_area_dir.glob(subj_pattern):
        if not p.is_dir():
            continue
        if re.fullmatch(r"sub-\d+", p.name):
            if p.name.split("sub-")[1].endswith(csv_id_digits):
                subject_candidates.append(p)

    if len(subject_candidates) == 0:
        raise FileNotFoundError(
            f"No subject directory found for CSV ID {csv_id_digits} "
            f"using pattern {subj_pattern} under {base_area_dir}"
        )
    if len(subject_candidates) > 1:
        raise ValueError(
            f"Multiple subject directories matched CSV ID {csv_id_digits}: "
            f"{sorted(p.name for p in subject_candidates)}"
        )

    subject_dir = subject_candidates[0]

    ses_pattern = f"ses-*{csv_ses_digits}"
    session_candidates = []

    for ses in subject_dir.glob(ses_pattern):
        if not ses.is_dir():
            continue
        if re.fullmatch(r"ses-\d+", ses.name):
            if ses.name.split("ses-")[1].endswith(csv_ses_digits):
                session_candidates.append(ses)

    if len(session_candidates) == 0:
        raise FileNotFoundError(
            f"No session directory found for CSV session ID {csv_ses_digits} "
            f"under {subject_dir.name} using pattern {ses_pattern}"
        )
    if len(session_candidates) > 1:
        raise ValueError(
            f"Multiple session directories matched CSV session ID {csv_ses_digits} "
            f"under {subject_dir.name}: {sorted(s.name for s in session_candidates)}"
        )

    return subject_dir.name, session_candidates[0].name


def find_subject_name_in_norm_dir(norm_dir: Path, csv_id_digits: str) -> str:
    """
    Resolve a PNC CSV subject ID directly against the normalized-loadings root.

    Example:
        CSV ID 992 -> <norm_dir>/sub-0992/
    """
    pattern = f"sub-*{csv_id_digits}"
    candidates = []

    for p in norm_dir.glob(pattern):
        if not p.is_dir():
            continue
        if re.fullmatch(r"sub-\d+", p.name):
            if p.name.split("sub-")[1].endswith(csv_id_digits):
                candidates.append(p)

    if len(candidates) == 0:
        raise FileNotFoundError(
            f"No normalized-loadings subject directory found for CSV ID "
            f"{csv_id_digits} under {norm_dir}"
        )
    if len(candidates) > 1:
        raise ValueError(
            f"Multiple normalized-loadings subject directories matched CSV ID "
            f"{csv_id_digits}: {sorted(p.name for p in candidates)}"
        )

    return candidates[0].name


def load_normed_csv(norm_csv_path: Path) -> np.ndarray:
    """
    Load a normalized PFN CSV and return float32 array shaped (59412, 17).

    The R script writes row.names = FALSE, so the CSV is expected to contain
    exactly 17 numeric PFN columns and no index column.
    """
    if not norm_csv_path.exists():
        raise FileNotFoundError(f"Normalized PFN CSV not found: {norm_csv_path}")

    df = pd.read_csv(norm_csv_path)

    if df.shape != (VERT_NUM, N_NETWORKS):
        raise ValueError(
            f"Unexpected normalized PFN CSV shape {df.shape} in {norm_csv_path}. "
            f"Expected ({VERT_NUM}, {N_NETWORKS})."
        )

    try:
        values = df.to_numpy(dtype=np.float32)
    except (TypeError, ValueError) as e:
        raise ValueError(
            f"Normalized PFN CSV contains non-numeric values: {norm_csv_path}"
        ) from e

    if not np.all(np.isfinite(values)):
        raise ValueError(
            f"Normalized PFN CSV contains non-finite values: {norm_csv_path}"
        )

    return values


def build_feature_cache(features_idx: np.ndarray):
    """
    Decode the PNC-derived flattened feature indices.

    Returns:
      - pos_by_net: positions in the final masked output vector for each PFN
      - vert_by_net: vertex indices corresponding to those positions
      - total_len: number of retained features
    """
    if features_idx.ndim != 1:
        raise ValueError(
            f"features_indices.npy must be 1D, got shape {features_idx.shape}"
        )

    if features_idx.size == 0:
        raise ValueError("features_indices.npy is empty")

    expected_max = N_NETWORKS * VERT_NUM
    if features_idx.min() < 0 or features_idx.max() >= expected_max:
        raise ValueError(
            f"Feature indices out of range. Expected [0, {expected_max - 1}], "
            f"got min={features_idx.min()} max={features_idx.max()}"
        )

    nets = (features_idx // VERT_NUM).astype(np.int16)
    vert = (features_idx % VERT_NUM).astype(np.int32)

    pos_by_net = []
    vert_by_net = []

    for network in range(N_NETWORKS):
        pos = np.where(nets == network)[0]
        pos_by_net.append(pos)
        vert_by_net.append(vert[pos])

    return pos_by_net, vert_by_net, features_idx.size


def extract_subject(
    sub_name: str,
    ses_name: str,
    pos_by_net,
    vert_by_net,
    total_len: int,
    features_idx: np.ndarray,
    net_id: np.ndarray,
    vertex_id: np.ndarray,
    par_out_dir: Path,
    norm_dir: Path,
    features_idx_path: Path,
) -> Path:
    """Read one normalized CSV, extract masked values, and write one .npz."""
    out_dir = par_out_dir / sub_name / ses_name / "PFN"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Keep the historical filename pattern so the existing compilation script
    # that searches for *PFN_loadings_masked_<N>.npz remains compatible.
    out_path = out_dir / f"{sub_name}_{ses_name}_PFN_loadings_masked_{total_len}.npz"

    if out_path.exists():
        print(f"{sub_name}_{ses_name}: output exists, skipping")
        return out_path

    norm_csv_path = norm_dir / sub_name / f"{sub_name}_PFN_loadings_normed.csv"
    pfn_loadings_norm = load_normed_csv(norm_csv_path)

    values = np.empty(total_len, dtype=np.float32)

    for network in range(N_NETWORKS):
        pos = pos_by_net[network]
        if pos.size == 0:
            continue

        v_idx = vert_by_net[network]
        values[pos] = pfn_loadings_norm[v_idx, network]

    np.savez_compressed(
        out_path,
        values=values,
        features_idx=features_idx.astype(np.int32),
        net_id=net_id,
        vertex_id=vertex_id,
        vert_num=np.int32(VERT_NUM),
        total_len=np.int32(total_len),
        indexing_note=np.array(
            "features_idx is 0-based and concatenated "
            "(network*vert_num + vertex0). net_id and vertex_id are 1-based "
            "for readability (PFN01..PFN17; vertex 1..vert_num).",
            dtype=object,
        ),
        normalization_note=np.array(
            "values are read from the precomputed R vertex-wise normalized PFN CSV; "
            "the 17 PFN loadings sum to 1 at each vertex when the original vertex "
            "sum was nonzero.",
            dtype=object,
        ),
        loading_source=np.array(str(norm_csv_path), dtype=object),
        features_idx_source=np.array(str(features_idx_path), dtype=object),
    )

    return out_path


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract PNC-mask-selected PFN loadings from precomputed normalized CSVs."
        )
    )

    parser.add_argument("--index", type=int, required=True)

    parser.add_argument(
        "--csv-path",
        type=str,
        default="/cbica/projects/bbl_22q/analysis/allometry/inputs/filtered_demo_22q.csv",
        help="Phenotype CSV. 22q requires subject/participant_id plus session; PNC requires subject/participant_id.",
    )

    parser.add_argument(
        "--norm-dir",
        type=str,
        default=None,
        help=(
            "Root containing <sub>/<sub>_PFN_loadings_normed.csv. "
            "Defaults to the 22q normalized directory, or the PNC normalized "
            "directory when --pnc is supplied."
        ),
    )

    parser.add_argument(
        "--base-area-dir",
        type=str,
        default="/cbica/projects/bbl_22q/analysis/allometry/outputs/22q_areas",
        help=(
            "22q area directory used only to resolve subject/session directory names. "
            "Not used for PNC subject resolution."
        ),
    )

    parser.add_argument(
        "--par-out-dir",
        type=str,
        default="/cbica/projects/bbl_22q/analysis/pfn/data/Final_PFNs_22q_nonzero_npz",
        help="Output root for masked normalized .npz files.",
    )

    parser.add_argument(
        "--feat-idx-path",
        type=str,
        default="/cbica/projects/bbl_22q/analysis/topography/inputs/features_indices.npy",
        help="Path to the PNC-derived nonzero feature indices .npy file.",
    )

    parser.add_argument(
        "--pnc",
        action="store_true",
        help="PNC mode: use ses-PNC1 and PNC normalized-loading defaults.",
    )

    args = parser.parse_args()

    csv_path = Path(args.csv_path)
    base_area_dir = Path(args.base_area_dir)
    par_out_dir = Path(args.par_out_dir)
    features_idx_path = Path(args.feat_idx_path)

    if args.norm_dir is not None:
        norm_dir = Path(args.norm_dir)
    elif args.pnc:
        norm_dir = Path(
            "/cbica/projects/bbl_22q/analysis/allometry/inputs/PNC_PFN_loadings_normed"
        )
    else:
        norm_dir = Path(
            "/cbica/projects/bbl_22q/analysis/allometry/inputs/22q_PFN_loadings_normed"
        )

    if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")
    if not norm_dir.exists():
        raise FileNotFoundError(f"Normalized PFN directory does not exist: {norm_dir}")
    if not args.pnc and not base_area_dir.exists():
        raise FileNotFoundError(f"Base area directory does not exist: {base_area_dir}")
    if not features_idx_path.exists():
        raise FileNotFoundError(
            f"features_indices.npy not found: {features_idx_path}"
        )

    par_out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(csv_path)

    if "subject" in df.columns:
        subj_col = "subject"
    elif "participant_id" in df.columns:
        subj_col = "participant_id"
    else:
        raise ValueError(
            f"CSV must contain a 'subject' or 'participant_id' column. "
            f"Found columns: {list(df.columns)}"
        )

    csv_ids = df[subj_col].tolist()

    if not args.pnc:
        if "session" not in df.columns:
            raise ValueError(
                f"22q CSV must contain a 'session' column. "
                f"Found columns: {list(df.columns)}"
            )
        csv_ses = df["session"].tolist()

    if len(csv_ids) == 0:
        raise ValueError(f"No rows found in CSV: {csv_path}")

    if args.index < 0 or args.index >= len(csv_ids):
        raise IndexError(
            f"--index out of range: {args.index} (n={len(csv_ids)})"
        )

    features_idx = np.load(features_idx_path, allow_pickle=False)
    pos_by_net, vert_by_net, total_len = build_feature_cache(features_idx)

    net_id = (features_idx // VERT_NUM).astype(np.int16) + 1
    vertex_id = (features_idx % VERT_NUM).astype(np.int32) + 1

    csv_id_digits = normalize_csv_id(csv_ids[args.index])

    if args.pnc:
        sub_name = find_subject_name_in_norm_dir(norm_dir, csv_id_digits)
        ses_name = "ses-PNC1"
        print(
            f"PNC mode: subject {csv_id_digits} -> {sub_name} / {ses_name}"
        )
    else:
        csv_ses_digits = normalize_csv_id(csv_ses[args.index])
        sub_name, ses_name = find_subject_dir_by_glob(
            base_area_dir,
            csv_id_digits,
            csv_ses_digits,
        )
        print(
            f"22q mode: subject {csv_id_digits}, session {csv_ses_digits} "
            f"-> {sub_name} / {ses_name}"
        )

    out_path = extract_subject(
        sub_name=sub_name,
        ses_name=ses_name,
        pos_by_net=pos_by_net,
        vert_by_net=vert_by_net,
        total_len=total_len,
        features_idx=features_idx,
        net_id=net_id,
        vertex_id=vertex_id,
        par_out_dir=par_out_dir,
        norm_dir=norm_dir,
        features_idx_path=features_idx_path,
    )

    print(f"Normalized loading source: {norm_dir / sub_name / f'{sub_name}_PFN_loadings_normed.csv'}")
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
