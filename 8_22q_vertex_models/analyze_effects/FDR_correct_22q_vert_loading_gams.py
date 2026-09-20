#!/usr/bin/env python3
"""
TCA-controlled:
python FDR_correct_22q_vert_loading_gams.py \
    --h5 "/cbica/projects/bbl_22q/analysis/topography/results/091426_normed_22q_vert_gam_results_nonstd_sex_avg_logTC_subclass.h5" \
    --outdir "/cbica/projects/bbl_22q/analysis/topography/results/091426_22q_logTC_vert_normed_loadings_FDR_pval_beta_subclass/"

no TCA cov:
python FDR_correct_22q_vert_loading_gams.py \
    --h5 "/cbica/projects/bbl_22q/analysis/topography/results/091426_normed_22q_vert_gam_results_nonstd_sex_avg_no_TCA_subclass.h5" \
    --outdir "/cbica/projects/bbl_22q/analysis/topography/results/091426_22q_no_TCA_vert_normed_loadings_FDR_pval_beta_subclass/"

Compute FDR corrected p-vals and save them out along with betas
"""

import argparse
import numpy as np
import pandas as pd
import h5py
from pathlib import Path
from scipy.stats import false_discovery_control

V = 59412
N_NET = 17
LH_N = 29696  # vertices 1..29696

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5", required=True, help="Path to compiled GAM results HDF5")
    parser.add_argument("--outdir", required=True, help="Output directory for csvs")
    args = parser.parse_args()

    h5_path = Path(args.h5)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with h5py.File(h5_path, "r") as f:

        # --- read in values ---
        net_id = f["net_id"][:].astype(np.int16)          # (TOTAL_LEN,)
        vertex_id = f["vertex_id"][:].astype(np.int32)    # (TOTAL_LEN,)
        beta = f["beta_stat_22q"][:].astype(np.float64)   # (TOTAL_LEN,)
        pval = f["p_stat_22q"][:].astype(np.float64)      # (TOTAL_LEN,)

        total_len = net_id.shape[0]
        if vertex_id.shape[0] != total_len or beta.shape[0] != total_len or pval.shape[0] != total_len:
            raise ValueError("Length mismatch among net_id/vertex_id/beta/pval")
        
        # Validate ranges
        if net_id.min() < 1 or net_id.max() > N_NET:
            raise ValueError(f"net_id out of range. min={net_id.min()} max={net_id.max()} expected [1,{N_NET}]")
        if vertex_id.min() < 1 or vertex_id.max() > V:
            raise ValueError(f"vertex_id out of range. min={vertex_id.min()} max={vertex_id.max()} expected [1,{V}]")

        # --- FDR correction of p-values, accounting for if any p values are NA ---
        fdr_p_vals = np.full_like(pval, np.nan)
        mask = np.isfinite(pval) & (pval >= 0) & (pval <= 1)
        print(f"Valid pvals for FDR: {mask.sum()} / {total_len}")
        if mask.any():
            fdr_p_vals[mask] = false_discovery_control(pval[mask], method="bh")

        # Convert vertex_id to 0-based indices for array assignment
        vert0 = vertex_id.astype(np.int64) - 1

        # ---- Write per-net CSVs: 59412 x 2 (beta, fdr_p), NaN for missing ----
        vertex_index = np.arange(1, V + 1, dtype=np.int32)  # 1..V for readability

        for net in range(1, N_NET + 1):
            net_features = (net_id == net)
            beta_net = np.full(V, np.nan, dtype=np.float64)
            fdr_net = np.full(V, np.nan, dtype=np.float64)

            if np.any(net_features):
                vert_idx = vert0[net_features]
                # Guard against any accidental out-of-range values
                in_range = (vert_idx >= 0) & (vert_idx < V)
                vert_idx = vert_idx[in_range]
                b_vals = beta[net_features][in_range]
                f_vals = fdr_p_vals[net_features][in_range]

                beta_net[vert_idx] = b_vals
                fdr_net[vert_idx] = f_vals

                # Guard against duplicate vertices
                if vert_idx.size != np.unique(vert_idx).size:
                    raise ValueError(f"Duplicate vertex entries for net {net}. Check compiled H5.")

            df = pd.DataFrame(
                {"beta": beta_net, "fdr_p": fdr_net},
                index=vertex_index
            )
            df.index.name = "vertex_id"

            out_path = outdir / f"stat_22q_gam_vert_loadings_PFN_{net:02d}_beta_fdr.csv"
            df.to_csv(out_path, na_rep="NA")
            print(f"Wrote {out_path}  (filled {np.isfinite(beta_net).sum()} / {V})")

    print("\nDone.")

if __name__ == "__main__":
    main()
