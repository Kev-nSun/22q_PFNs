#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
from typing import Dict, List, Optional

import h5py
import numpy as np
import pandas as pd
import scipy.io
import os
from scipy.io import savemat
from scipy.stats import median_abs_deviation


"""
Calculate the median absolute deviation of vertex-wise normalized PFN
loadings across QC-passed participants (based on --qc-csv provided).

For each participant and each cortical vertex, the 17 PFN loadings are
normalized to sum to 1:

    normalized_loading[v, k] =
        loading[v, k] / sum(loading[v, 0:17])

Any NaN or infinite values produced by normalization are replaced with zero,
matching the associated R normalization script.

MAD is then calculated across participants separately for every one of the
1,010,004 network-vertex combinations.

No nonzero feature mask is used.
"""


VERT_NUM = 59412
NETWORK_NUM = 17
TOTAL_FEATURES = VERT_NUM * NETWORK_NUM  # 1,010,004

run_id = os.environ.get("SLURM_JOB_ID", str(os.getpid()))

def normalize_csv_id(value) -> str:
    """
    Normalize a subject identifier to digits.

    Examples:
        sub-0992 -> 0992
        992.0    -> 992
    """
    subject_id = str(value).strip()
    subject_id = subject_id.replace("sub-", "")
    subject_id = re.sub(r"\.0$", "", subject_id)

    if not re.fullmatch(r"\d+", subject_id):
        raise ValueError(
            f"Subject ID is not purely numeric after cleaning: "
            f"{value!r} -> {subject_id!r}"
        )

    return subject_id


def choose_fn_array(arrays: Dict[str, np.ndarray], source: Path) -> np.ndarray:
    """
    Select the PFN matrix from arrays read from a MATLAB file.

    Preference:
      1. A variable named FN
      2. A unique numeric 2-D array with one accepted FN shape
    """
    if "FN" in arrays:
        candidate = np.asarray(arrays["FN"])

        if candidate.ndim == 2 and candidate.shape in {
            (VERT_NUM, NETWORK_NUM),
            (NETWORK_NUM, VERT_NUM),
        }:
            return candidate

        raise ValueError(
            f"Variable 'FN' in {source} has shape {candidate.shape}; "
            f"expected ({VERT_NUM}, {NETWORK_NUM}) or "
            f"({NETWORK_NUM}, {VERT_NUM})."
        )

    candidates = []

    for key, value in arrays.items():
        if key.startswith("__"):
            continue

        value = np.asarray(value)

        if (
            value.ndim == 2
            and np.issubdtype(value.dtype, np.number)
            and value.shape
            in {
                (VERT_NUM, NETWORK_NUM),
                (NETWORK_NUM, VERT_NUM),
            }
        ):
            candidates.append((key, value))

    if len(candidates) == 1:
        return candidates[0][1]

    if len(candidates) > 1:
        candidate_names = [key for key, _ in candidates]
        raise ValueError(
            f"Multiple plausible PFN matrices found in {source}: "
            f"{candidate_names}. Rename the desired variable to 'FN'."
        )

    raise ValueError(
        f"No valid PFN matrix was found in {source}. Expected a numeric "
        f"array shaped ({VERT_NUM}, {NETWORK_NUM}) or "
        f"({NETWORK_NUM}, {VERT_NUM})."
    )


def load_v73_mat(path: Path) -> np.ndarray:
    """
    Read a MATLAB v7.3 HDF5 file.

    MATLAB/HDF5 dimension ordering can differ from scipy.io.loadmat output,
    so both accepted orientations are checked after loading.
    """
    arrays = {}

    with h5py.File(path, "r") as mat_file:
        if "FN" in mat_file and isinstance(mat_file["FN"], h5py.Dataset):
            arrays["FN"] = np.asarray(mat_file["FN"])
        else:
            for key, value in mat_file.items():
                if isinstance(value, h5py.Dataset):
                    arrays[key] = np.asarray(value)

    return choose_fn_array(arrays, path)


def load_fn_mat(path: Path) -> np.ndarray:
    """
    Load an FN.mat file and return an array shaped (59412, 17).

    Supports traditional MATLAB files with scipy.io.loadmat and MATLAB
    v7.3 HDF5 files with h5py.
    """
    if not path.exists():
        raise FileNotFoundError(f"FN.mat not found: {path}")

    try:
        mat_data = scipy.io.loadmat(str(path))
        fn = choose_fn_array(mat_data, path)
    except (NotImplementedError, ValueError) as scipy_error:
        # Try HDF5 only when the file is actually an HDF5 file.
        if not h5py.is_hdf5(path):
            raise scipy_error

        fn = load_v73_mat(path)

    fn = np.asarray(fn)

    if fn.shape == (VERT_NUM, NETWORK_NUM):
        normalized = fn
    elif fn.shape == (NETWORK_NUM, VERT_NUM):
        normalized = fn.T
    else:
        raise ValueError(
            f"Unexpected FN shape {fn.shape} in {path}. Expected "
            f"({VERT_NUM}, {NETWORK_NUM}) or "
            f"({NETWORK_NUM}, {VERT_NUM})."
        )

    if not np.issubdtype(normalized.dtype, np.number):
        raise TypeError(f"FN matrix is not numeric in {path}")

    return normalized


def find_subject_directory(fn_dir: Path, csv_id_digits: str) -> Path:
    """
    Match one CSV subject ID to one sub-* directory by numeric suffix.

    This preserves the behavior used by the earlier extraction script. For
    example, CSV ID 992 may match sub-0992.
    """
    pattern = f"sub-*{csv_id_digits}"

    candidates = [
        path
        for path in fn_dir.glob(pattern)
        if path.is_dir()
        and re.fullmatch(r"sub-\d+", path.name)
        and path.name.removeprefix("sub-").endswith(csv_id_digits)
    ]

    if len(candidates) == 0:
        raise FileNotFoundError(
            f"No subject directory matched CSV ID {csv_id_digits} "
            f"under {fn_dir}"
        )

    if len(candidates) > 1:
        raise ValueError(
            f"Multiple subject directories matched CSV ID {csv_id_digits}: "
            f"{sorted(path.name for path in candidates)}"
        )

    return candidates[0]


def get_fn_paths(
    fn_dir: Path,
    fn_filename: str,
    qc_csv: Optional[Path],
) -> List[Path]:
    """
    Obtain participant FN.mat paths.

    If a CSV is supplied, its subject order and participant selection are used.
    Otherwise, all matching sub-*/FN.mat files are used.
    """
    if qc_csv is None:
        paths = sorted(fn_dir.glob(f"sub-*/{fn_filename}"))

        if not paths:
            raise FileNotFoundError(
                f"No files matched {fn_dir}/sub-*/{fn_filename}"
            )

        return paths

    if not qc_csv.exists():
        raise FileNotFoundError(f"CSV file not found: {qc_csv}")

    dataframe = pd.read_csv(qc_csv, dtype=str)
    dataframe.columns = dataframe.columns.str.strip()

    if "subject" in dataframe.columns:
        subject_column = "subject"
    elif "participant_id" in dataframe.columns:
        subject_column = "participant_id"
    else:
        raise ValueError(
            "CSV must contain a 'subject' or 'participant_id' column. "
            f"Found: {list(dataframe.columns)}"
        )

    fn_paths = []
    seen_subjects = set()

    for row_number, raw_subject_id in enumerate(
        dataframe[subject_column].tolist()
    ):
        subject_digits = normalize_csv_id(raw_subject_id)
        subject_dir = find_subject_directory(fn_dir, subject_digits)

        # FN.mat has no session component in the specified directory structure.
        # Avoid processing duplicate CSV rows for the same participant.
        if subject_dir.name in seen_subjects:
            print(
                f"Row {row_number}: duplicate subject {subject_dir.name}; "
                "skipping duplicate."
            )
            continue

        fn_path = subject_dir / fn_filename

        if not fn_path.exists():
            raise FileNotFoundError(
                f"FN file missing for {subject_dir.name}: {fn_path}"
            )

        seen_subjects.add(subject_dir.name)
        fn_paths.append(fn_path)

    if not fn_paths:
        raise ValueError(f"No participants were resolved from {qc_csv}")

    return fn_paths


def normalize_fn_per_vertex(fn: np.ndarray) -> np.ndarray:
    """
    Normalize each vertex's 17 PFN loadings so they sum to 1.

    This reproduces the R operations:

        vertex_sums <- rowSums(pfn_loadings)
        pfn_loadings_norm <- pfn_loadings / vertex_sums
        pfn_loadings_norm[!is.finite(pfn_loadings_norm)] <- 0

    Parameters
    ----------
    fn
        PFN matrix shaped (59412 vertices, 17 networks).

    Returns
    -------
    np.ndarray
        Float32 matrix shaped (59412, 17), normalized across networks
        separately for each vertex.
    """
    if fn.shape != (VERT_NUM, NETWORK_NUM):
        raise ValueError(
            f"Internal FN shape error: got {fn.shape}, expected "
            f"({VERT_NUM}, {NETWORK_NUM})"
        )

    fn = np.asarray(fn, dtype=np.float32)

    # One sum per cortical vertex across the 17 networks.
    vertex_sums = np.sum(fn, axis=1, keepdims=True)

    # Match the R behavior: division by zero may produce NaN or Inf,
    # which is subsequently replaced with zero.
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        fn_norm = fn / vertex_sums

    fn_norm[~np.isfinite(fn_norm)] = 0.0

    return fn_norm.astype(np.float32, copy=False)


def flatten_fn_network_major(fn: np.ndarray) -> np.ndarray:
    """
    Normalize PFN loadings across the 17 networks at each vertex, then
    flatten in network-major order.

    Output order:

        PFN01 vertices 1..59412,
        PFN02 vertices 1..59412,
        ...
        PFN17 vertices 1..59412

    The flattened index is:

        network_zero_based * 59412 + vertex_zero_based
    """
    fn_norm = normalize_fn_per_vertex(fn)

    return fn_norm.T.reshape(TOTAL_FEATURES)


def create_feature_memmap(
    fn_paths: List[Path],
    memmap_path: Path,
    storage_dtype: np.dtype,
) -> np.memmap:
    """
    Read every participant's FN.mat, normalize each vertex across the
    17 PFNs, and create a participant-by-feature disk-backed matrix.
    """
    participant_count = len(fn_paths)

    feature_matrix = np.memmap(
        memmap_path,
        mode="w+",
        dtype=storage_dtype,
        shape=(participant_count, TOTAL_FEATURES),
    )

    for participant_index, fn_path in enumerate(fn_paths):
        subject_name = fn_path.parent.name

        print(
            f"[{participant_index + 1}/{participant_count}] "
            f"Loading and normalizing {subject_name}: {fn_path}"
        )

        fn = load_fn_mat(fn_path)

        # Normalization now occurs inside this function.
        flattened_normed = flatten_fn_network_major(fn)

        if not np.all(np.isfinite(flattened_normed)):
            raise RuntimeError(
                f"Non-finite normalized values remain for {fn_path}. "
                "This should not occur because non-finite values are "
                "replaced with zero during normalization."
            )

        feature_matrix[participant_index, :] = (
            flattened_normed.astype(storage_dtype, copy=False)
        )

    feature_matrix.flush()
    return feature_matrix


def calculate_mad(
    feature_matrix: np.memmap,
    chunk_size: int,
) -> np.ndarray:
    """
    Calculate exact median absolute deviation along the participant axis.

    scipy.stats.median_abs_deviation uses:

        median(abs(x - median(x)))

    The default scale=1 is retained, so this is raw MAD rather than a
    normal-distribution-scaled estimate.
    """
    if chunk_size <= 0:
        raise ValueError("--chunk-size must be greater than zero")

    mad_values = np.empty(TOTAL_FEATURES, dtype=np.float32)

    for start in range(0, TOTAL_FEATURES, chunk_size):
        stop = min(start + chunk_size, TOTAL_FEATURES)

        print(
            f"Calculating MAD for features "
            f"{start:,} through {stop - 1:,} "
            f"of {TOTAL_FEATURES:,}"
        )

        # Convert the current chunk to float32 for calculation. This keeps
        # memory bounded while avoiding float16 arithmetic during the MAD.
        chunk = np.asarray(
            feature_matrix[:, start:stop],
            dtype=np.float32,
        )

        chunk_mad = median_abs_deviation(
            chunk,
            axis=0,
            scale=1.0,
            nan_policy="propagate",
        )

        mad_values[start:stop] = chunk_mad.astype(np.float32, copy=False)

    return mad_values


def save_results(
    output_dir: Path,
    mad_values: np.ndarray,
    participant_count: int,
    save_npy: bool,
) -> Path:
    """
    Save MAD values calculated from vertex-wise normalized PFN loadings.
    """
    output_path = output_dir / "PFN_loading_normed_MAD.mat"

    flat_indices = np.arange(TOTAL_FEATURES, dtype=np.int32)

    net_id = (flat_indices // VERT_NUM) + 1
    vertex_id = (flat_indices % VERT_NUM) + 1

    savemat(
        output_path,
        {
            "MAD": mad_values.reshape(1, -1),
            "net_id": net_id.reshape(1, -1),
            "vertex_id": vertex_id.reshape(1, -1),
            "vert_num": np.array([[VERT_NUM]], dtype=np.int32),
            "network_num": np.array([[NETWORK_NUM]], dtype=np.int32),
            "total_features": np.array(
                [[TOTAL_FEATURES]],
                dtype=np.int32,
            ),
            "participant_count": np.array(
                [[participant_count]],
                dtype=np.int32,
            ),
            "normalization_applied": np.array(
                [["vertex-wise sum-to-one across the 17 PFNs"]],
                dtype=object,
            ),
            "normalization_definition": np.array(
                [
                    [
                        "For each participant and cortical vertex: "
                        "FN_norm[vertex, network] = "
                        "FN[vertex, network] / sum(FN[vertex, 1:17]). "
                        "NaN and Inf values were replaced with zero."
                    ]
                ],
                dtype=object,
            ),
            "indexing_note": np.array(
                [
                    [
                        "MAD is stored in network-major order: "
                        "PFN01 vertices 1..59412, followed by PFN02, "
                        "through PFN17. No nonzero feature mask was used."
                    ]
                ],
                dtype=object,
            ),
            "mad_definition": np.array(
                [
                    [
                        "median(abs(x - median(x))) across "
                        "participants, calculated from vertex-wise "
                        "normalized PFN loadings; scipy scale=1.0"
                    ]
                ],
                dtype=object,
            ),
        },
        do_compression=True,
    )

    if save_npy:
        np.save(
            output_dir / "PFN_loading_normed_MAD.npy",
            mad_values,
            allow_pickle=False,
        )

    return output_path

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate MAD across participants for all 1,010,004 "
            "vertex-wise normalized PFN network-vertex loadings. "
            "No nonzero feature mask is used."
        )
    )

    parser.add_argument(
        "output_dir",
        type=Path,
        help="Directory in which results and temporary files are written.",
    )
    parser.add_argument(
        "fn_dir",
        type=Path,
        help="Directory containing sub-*/FN.mat.",
    )
    parser.add_argument(
        "--qc-csv",
        type=Path,
        default=None,
        required=True,
        help=(
            "Optional participant CSV containing 'subject' or "
            "'participant_id'. If omitted, all sub-*/FN.mat files are used."
        ),
    )
    parser.add_argument(
        "--fn-filename",
        default="FN.mat",
        help="FN filename within each subject directory. Default: FN.mat",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=4096,
        help=(
            "Number of features processed per MAD chunk. Lower this if "
            "memory is limited. Default: 4096"
        ),
    )
    parser.add_argument(
        "--temp-dir",
        type=Path,
        default=None,
        help=(
            "Directory for the temporary participant-by-feature matrix. "
            "Defaults to output_dir."
        ),
    )
    parser.add_argument(
        "--keep-memmap",
        action="store_true",
        help="Keep the temporary participant-by-feature memmap.",
    )
    parser.add_argument(
        "--save-npy",
        action="store_true",
        help="Also save PFN_loading_normed_MAD.npy.",
    )

    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    fn_dir = args.fn_dir.resolve()
    qc_csv = args.qc_csv.resolve()

    if not fn_dir.exists():
        raise FileNotFoundError(f"FN directory not found: {fn_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    fn_paths = get_fn_paths(
        fn_dir=fn_dir,
        fn_filename=args.fn_filename,
        qc_csv=qc_csv,
    )

    participant_count = len(fn_paths)
    storage_dtype = np.dtype("float32")
    temp_dir = (
        args.temp_dir.resolve()
        if args.temp_dir is not None
        else output_dir
    )
    temp_dir.mkdir(parents=True, exist_ok=True)

    memmap_path = temp_dir / (
        f"PFN_all_participants_"
        f"{participant_count}x{TOTAL_FEATURES}_"
        f"{storage_dtype.name}_{run_id}.dat"
    )

    bytes_required = (
        participant_count
        * TOTAL_FEATURES
        * storage_dtype.itemsize
    )

    print(f"Participants: {participant_count:,}")
    print(f"Features per participant: {TOTAL_FEATURES:,}")
    print(f"Temporary matrix: {memmap_path}")
    print(
        f"Approximate temporary disk requirement: "
        f"{bytes_required / (1024 ** 3):.2f} GiB"
    )
    print("No nonzero feature mask will be applied.")

    feature_matrix = None

    try:
        feature_matrix = create_feature_memmap(
            fn_paths=fn_paths,
            memmap_path=memmap_path,
            storage_dtype=storage_dtype,
        )

        mad_values = calculate_mad(
            feature_matrix=feature_matrix,
            chunk_size=args.chunk_size,
        )

        output_path = save_results(
            output_dir=output_dir,
            mad_values=mad_values,
            participant_count=participant_count,
            save_npy=args.save_npy,
        )

        print(f"Saved MAD results: {output_path}")

    finally:
        # Release the memory map before attempting deletion.
        if feature_matrix is not None:
            feature_matrix.flush()
            del feature_matrix

        if memmap_path.exists() and not args.keep_memmap:
            memmap_path.unlink()
            print(f"Removed temporary matrix: {memmap_path}")


if __name__ == "__main__":
    main()