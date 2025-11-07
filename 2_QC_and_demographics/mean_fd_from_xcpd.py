# python mean_fd_from_xcpd.py --base /cbica/projects/bbl_22q/data/22q_xcpd_extracted --out /cbica/projects/bbl_22q/data/22q_xcpd_extracted/22q_mean_fd_summary.csv
# python mean_fd_from_xcpd.py --base /cbica/projects/bbl_22q/data/PNC_xcpd_extracted/xcpd --out /cbica/projects/bbl_22q/data/PNC_xcpd_extracted/PNC_mean_fd_summary.csv

#!/usr/bin/env python3
import os, glob, argparse
import pandas as pd
import numpy as np

# Read just FD column from one motion file
def read_fd(fp):
    if not fp or not os.path.isfile(fp):
        return pd.Series(dtype=float)
    try:
        df = pd.read_csv(fp, sep="\t", usecols=["framewise_displacement"])
        s = pd.to_numeric(df["framewise_displacement"], errors="coerce")
    except ValueError:
        return pd.Series(dtype=float)
    return s.replace([np.inf, -np.inf], np.nan).dropna()
    
def main():
    parser = argparse.ArgumentParser(
        description="Compute mean FD from XCP-D motion.tsv files for idemo/rest and combined."
    )
    parser.add_argument(
        "--base",
        required=True,
        help="Dataset directory containing sub-*/ses-*/func/*motion.tsv",
    )
    parser.add_argument(
        "-o", "--out",
        default=None,
        help="Optional output CSV path (if omitted, prints to stdout)."
    )
    args = parser.parse_args()

    BASE = args.base

    # Find motion files for both tasks
    idemo_paths = glob.glob(os.path.join(BASE, "sub-*", "ses-*", "func", "*task-idemo_*motion.tsv"))
    rest_paths  = glob.glob(os.path.join(BASE, "sub-*", "ses-*", "func", "*task-rest_*motion.tsv"))

    def key_from_path(p):
        # expects .../sub-XXX/ses-YYY/func/...motion.tsv
        parts = p.split(os.sep)
        sub = next(x for x in parts if x.startswith("sub-"))
        ses = next(x for x in parts if x.startswith("ses-"))
        return sub, ses

    # Map (sub, ses) -> file path per task
    records = {}
    for p in idemo_paths:
        records.setdefault(key_from_path(p), {})["idemo"] = p
    for p in rest_paths:
        records.setdefault(key_from_path(p), {})["rest"] = p

    rows = []
    for (sub, ses), tasks in sorted(records.items()):
        idemo_path = tasks.get("idemo", "")
        rest_path  = tasks.get("rest", "")

        fd_idemo = read_fd(idemo_path) if idemo_path else pd.Series(dtype=float)
        fd_rest  = read_fd(rest_path)  if rest_path  else pd.Series(dtype=float)

        mean_idemo = fd_idemo.mean() if not fd_idemo.empty else float("nan")
        mean_rest  = fd_rest.mean()  if not fd_rest.empty  else float("nan")

        # Combined = mean FD across concatenated volumes of both scans (or whatever exists)
        parts = [s for s in (fd_idemo, fd_rest) if not s.empty]
        fd_all = pd.concat(parts, ignore_index=True) if parts else pd.Series(dtype=float)
        mean_combined = fd_all.mean() if not fd_all.empty else float("nan")

        rows.append({
            "sub": sub,
            "ses": ses,
            "n_vols_idemo": int(fd_idemo.size),
            "mean_fd_idemo": mean_idemo,
            "n_vols_rest": int(fd_rest.size),
            "mean_fd_rest": mean_rest,
            "n_vols_total": int(fd_all.size),
            "mean_fd_combined": mean_combined,
            "path_idemo": tasks.get("idemo", ""),
            "path_rest": tasks.get("rest", ""),
        })

    out = pd.DataFrame(rows).sort_values(["sub", "ses"]).reset_index(drop=True)

    if args.out:
        out.to_csv(args.out, index=False)
        print(f"Wrote: {args.out} ({len(out)} rows)")
    else:
        print(out.to_csv(index=False))

if __name__ == "__main__":
    main()