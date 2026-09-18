#!/usr/bin/env python3
"""
python hard_parcellate_loading_gam_predictions_100_frames_normed.py \
  --h5 "/cbica/projects/bbl_22q/analysis/topography/results/091426_normed_TCA_vert_gam_results_nonstd_sex_avg_logTC_compiled.h5" \
  --outdir "/cbica/projects/bbl_22q/analysis/topography/results/091426_normed_gam_hard_parcellation_frames"

Compute winner (hard) PFN parcellation per cortical vertex for a range of TC
area values based on compiled GAM predictions stored in an HDF5.

Builds the mapping using:
  - net_id    (1..17)
  - vertex_id (1..59412)

Non-finite GAM predictions (NaN, +Inf, -Inf) are treated as ineligible
for winning the hard parcellation.

Outputs:
  Frame_<###>_PFN_GAM_hard_parcellation_L.mat
  Frame_<###>_PFN_GAM_hard_parcellation_R.mat

Each MAT file contains:
  hard_parcellation : int16 vector (n_vertices x 1), values 0..17
                      (0 means no eligible networks)
"""

import argparse
from pathlib import Path

import h5py
import numpy as np
from scipy.io import savemat


V = 59412
N_NET = 17
LH_N = 29696
RH_N = V - LH_N  # 29716


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--h5",
        required=True,
        help="Path to compiled GAM results HDF5",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for MAT frames",
    )
    parser.add_argument(
        "--block",
        type=int,
        default=4096,
        help="Vertex block size (default: 4096)",
    )
    args = parser.parse_args()

    h5_path = Path(args.h5)
    outdir = Path(args.outdir)

    if args.block <= 0:
        raise ValueError(f"--block must be > 0, got {args.block}")

    if not h5_path.exists():
        raise FileNotFoundError(f"HDF5 file not found: {h5_path}")

    outdir.mkdir(parents=True, exist_ok=True)

    with h5py.File(h5_path, "r") as f:

        # ------------------------------------------------------------
        # Validate required identifier datasets
        # ------------------------------------------------------------
        if "net_id" not in f:
            raise KeyError("net_id not found in HDF5")

        if "vertex_id" not in f:
            raise KeyError("vertex_id not found in HDF5")

        # net_id and vertex_id are 1-based in the HDF5
        net_id = f["net_id"][:].astype(np.int16)
        vertex_id = f["vertex_id"][:].astype(np.int32)

        total_len = net_id.shape[0]

        if vertex_id.shape[0] != total_len:
            raise ValueError(
                "net_id and vertex_id lengths differ: "
                f"{total_len} vs {vertex_id.shape[0]}"
            )

        if total_len == 0:
            raise ValueError("net_id / vertex_id datasets are empty")

        # ------------------------------------------------------------
        # Validate identifier ranges
        # ------------------------------------------------------------
        if net_id.min() < 1 or net_id.max() > N_NET:
            raise ValueError(
                f"net_id out of range: "
                f"min={net_id.min()} max={net_id.max()} "
                f"expected [1,{N_NET}]"
            )

        if vertex_id.min() < 1 or vertex_id.max() > V:
            raise ValueError(
                f"vertex_id out of range: "
                f"min={vertex_id.min()} max={vertex_id.max()} "
                f"expected [1,{V}]"
            )

        print(f"TOTAL_LEN = {total_len}")

        # Convert IDs to 0-based indices.
        net0 = net_id.astype(np.int32) - 1
        vert0 = vertex_id.astype(np.int64) - 1

        # ------------------------------------------------------------
        # Validate uniqueness of (vertex_id, net_id)
        # ------------------------------------------------------------
        #
        # Each vertex/network pair must correspond to at most one row.
        # Without this check, duplicate pairs would be silently
        # overwritten when rows_by_vertex is constructed.
        #
        flat_key = vert0 * N_NET + net0
        n_unique_pairs = np.unique(flat_key).size

        if n_unique_pairs != total_len:
            n_duplicate_rows = total_len - n_unique_pairs
            raise ValueError(
                "Duplicate (vertex_id, net_id) entries detected: "
                f"{n_duplicate_rows} duplicate row(s)"
            )

        print("Validated unique (vertex_id, net_id) pairs.")

        # ------------------------------------------------------------
        # Locate prediction dataset
        # ------------------------------------------------------------
        if "pred" not in f or "pred" not in f["pred"]:
            if "pred" in f:
                found = list(f["pred"].keys())
                raise KeyError(
                    f"pred/pred not found in HDF5. "
                    f"Found under pred/: {found}"
                )
            raise KeyError("pred group not found in HDF5")

        pred_ds = f["pred"]["pred"]

        if pred_ds.ndim != 2:
            raise ValueError(
                f"pred/pred must be 2-D, got shape {pred_ds.shape}"
            )

        pred_shape = pred_ds.shape

        # ------------------------------------------------------------
        # Determine grid_n robustly
        # ------------------------------------------------------------
        if "meta" in f and "grid_n" in f["meta"]:
            grid_n = int(np.asarray(f["meta"]["grid_n"]).ravel()[0])

            # Determine orientation using declared grid_n.
            if pred_shape == (total_len, grid_n):
                pred_mode = "row_major"
            elif pred_shape == (grid_n, total_len):
                pred_mode = "frame_major"
            else:
                raise ValueError(
                    f"Unexpected pred/pred shape {pred_shape}. "
                    f"With TOTAL_LEN={total_len} and grid_n={grid_n}, "
                    f"expected ({total_len}, {grid_n}) or "
                    f"({grid_n}, {total_len})"
                )

        else:
            # Infer the frame axis from TOTAL_LEN.
            if pred_shape[0] == total_len and pred_shape[1] != total_len:
                grid_n = pred_shape[1]
                pred_mode = "row_major"

            elif pred_shape[1] == total_len and pred_shape[0] != total_len:
                grid_n = pred_shape[0]
                pred_mode = "frame_major"

            elif pred_shape[0] == total_len and pred_shape[1] == total_len:
                raise ValueError(
                    "Cannot infer prediction orientation because both "
                    f"dimensions equal TOTAL_LEN={total_len}. "
                    "Provide meta/grid_n in the HDF5."
                )

            else:
                raise ValueError(
                    f"Cannot infer grid_n or prediction orientation from "
                    f"pred/pred shape {pred_shape} with "
                    f"TOTAL_LEN={total_len}"
                )

        if grid_n <= 0:
            raise ValueError(f"Invalid grid_n={grid_n}")

        print(f"grid_n = {grid_n}")
        print(
            f"pred/pred shape = {pred_shape}, "
            f"dtype = {pred_ds.dtype}, "
            f"using mode: {pred_mode}"
        )

        # ------------------------------------------------------------
        # Build (vertex, network) -> HDF5 row mapping
        # ------------------------------------------------------------
        #
        # rows_by_vertex[v0, net0] =
        #     masked row index (0..TOTAL_LEN-1)
        # or -1 if that network is unavailable for the vertex.
        #
        print(
            "Building (vertex, net) -> masked row mapping "
            "from net_id + vertex_id ..."
        )

        rows_by_vertex = np.full(
            (V, N_NET),
            -1,
            dtype=np.int32,
        )

        rows = np.arange(total_len, dtype=np.int32)
        rows_by_vertex[vert0, net0] = rows

        n_available_pairs = np.count_nonzero(rows_by_vertex >= 0)
        n_missing_pairs = rows_by_vertex.size - n_available_pairs

        print(
            f"Available vertex/network pairs: {n_available_pairs}"
        )
        print(
            f"Missing vertex/network pairs:   {n_missing_pairs}"
        )

        # ------------------------------------------------------------
        # winner[v0, frame] =
        #     winning network ID 1..17
        #     0 if no network has a finite prediction
        # ------------------------------------------------------------
        winner = np.zeros((V, grid_n), dtype=np.int16)

        block = int(args.block)

        # Aggregate prediction QC counters.
        total_nan = 0
        total_posinf = 0
        total_neginf = 0

        print("Computing hard parcellation...")

        for v0_start in range(0, V, block):
            v0_end = min(V, v0_start + block)
            vb = v0_end - v0_start

            rows_block = rows_by_vertex[v0_start:v0_end, :]
            # shape: (vb, 17)

            # Use float64 so argmax is based on the highest practical
            # precision available from typical GAM prediction datasets.
            #
            # Missing/ineligible networks remain -Inf.
            pred_block = np.full(
                (vb, N_NET, grid_n),
                -np.inf,
                dtype=np.float64,
            )

            for k in range(N_NET):
                rows_k = rows_block[:, k]
                valid = rows_k >= 0

                if not np.any(valid):
                    continue

                if pred_mode == "row_major":
                    vals = pred_ds[rows_k[valid], :]
                else:
                    vals = pred_ds[:, rows_k[valid]].T

                vals = np.asarray(vals, dtype=np.float64)

                # --------------------------------------------
                # Diagnose non-finite GAM predictions
                # --------------------------------------------
                n_nan = int(np.count_nonzero(np.isnan(vals)))
                n_posinf = int(np.count_nonzero(np.isposinf(vals)))
                n_neginf = int(np.count_nonzero(np.isneginf(vals)))

                total_nan += n_nan
                total_posinf += n_posinf
                total_neginf += n_neginf

                if n_nan or n_posinf or n_neginf:
                    print(
                        f"WARNING vertices "
                        f"{v0_start + 1}..{v0_end}, "
                        f"net {k + 1}: "
                        f"NaN={n_nan}, "
                        f"+Inf={n_posinf}, "
                        f"-Inf={n_neginf}"
                    )

                # --------------------------------------------
                # Any non-finite prediction is ineligible.
                # --------------------------------------------
                vals[~np.isfinite(vals)] = -np.inf

                pred_block[valid, k, :] = vals

            # ------------------------------------------------
            # Hard network winner across the 17 networks
            # ------------------------------------------------
            max_pred = np.max(pred_block, axis=1)
            all_ineligible = np.isneginf(max_pred)

            argmax_k = (
                np.argmax(pred_block, axis=1).astype(np.int16) + 1
            )

            # Label 0 means no network had a finite prediction.
            argmax_k[all_ineligible] = 0

            winner[v0_start:v0_end, :] = argmax_k

            if (
                (v0_start // block) % 5 == 0
                or v0_end == V
            ):
                print(
                    f"  vertices "
                    f"{v0_start + 1}..{v0_end} / {V}"
                )

        # ------------------------------------------------------------
        # Aggregate QC
        # ------------------------------------------------------------
        print("\nPrediction QC:")
        print(f"  NaN predictions:  {total_nan}")
        print(f"  +Inf predictions: {total_posinf}")
        print(f"  -Inf predictions: {total_neginf}")

        total_nonfinite = (
            total_nan + total_posinf + total_neginf
        )

        print(
            f"  Total non-finite:  {total_nonfinite}"
        )

        n_zero = int(np.count_nonzero(winner == 0))
        n_winner_values = winner.size

        zero_pct = (
            100.0 * n_zero / n_winner_values
            if n_winner_values
            else 0.0
        )

        print("\nHard-parcellation QC:")
        print(
            f"  Label-0 vertex/frame positions: "
            f"{n_zero} / {n_winner_values} "
            f"({zero_pct:.6f}%)"
        )

        # Sanity-check output label range.
        winner_min = int(winner.min())
        winner_max = int(winner.max())

        print(
            f"  Winner label range: "
            f"{winner_min}..{winner_max}"
        )

        if winner_min < 0 or winner_max > N_NET:
            raise RuntimeError(
                f"Unexpected winner labels: "
                f"min={winner_min}, max={winner_max}"
            )

        # ------------------------------------------------------------
        # Save per-frame MAT files for LH and RH
        # ------------------------------------------------------------
        print("\nSaving MAT frames...")

        for frame in range(grid_n):
            frame_num = frame + 1
            frame_str = f"{frame_num:03d}"

            lh_vec = winner[:LH_N, frame].reshape(-1, 1)
            rh_vec = winner[LH_N:, frame].reshape(-1, 1)

            # Explicit shape sanity checks.
            if lh_vec.shape != (LH_N, 1):
                raise RuntimeError(
                    f"Unexpected LH shape: {lh_vec.shape}"
                )

            if rh_vec.shape != (RH_N, 1):
                raise RuntimeError(
                    f"Unexpected RH shape: {rh_vec.shape}"
                )

            lh_name = (
                f"Frame_{frame_str}_"
                f"PFN_GAM_hard_parcellation_L.mat"
            )

            rh_name = (
                f"Frame_{frame_str}_"
                f"PFN_GAM_hard_parcellation_R.mat"
            )

            savemat(
                outdir / lh_name,
                {"hard_parcellation": lh_vec},
                do_compression=True,
            )

            savemat(
                outdir / rh_name,
                {"hard_parcellation": rh_vec},
                do_compression=True,
            )

            if frame_num % 10 == 0 or frame_num == grid_n:
                print(
                    f"  wrote frames "
                    f"{frame_num}/{grid_n}"
                )

    print("\nDone.")


if __name__ == "__main__":
    main()