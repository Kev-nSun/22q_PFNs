#!/usr/bin/env python3
"""
merge_dtseries_by_session.py

Concatenate dtseries within each (sub, ses) using Connectome Workbench
(wb_command -cifti-merge) and compile a global scan list of merged outputs for pNet.

Example:
  python merge_dtseries_by_session.py \
    --root /cbica/projects/PFN_ABCD/abcd-hcp-pipeline_0.1.4_timeseries \
    --out  /cbica/projects/PFN_ABCD/long_PFN_scripts/pnet_inputs/merged_dtseries
"""

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

BIDS_SUB_RE = re.compile(r"(sub-[A-Za-z0-9]+)")
BIDS_SES_RE = re.compile(r"(ses-[A-Za-z0-9]+)")

def which(cmd: str) -> str:
    p = shutil.which(cmd)
    if not p:
        sys.exit(f"[ERROR] '{cmd}' not found in PATH.")
    return p

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Concatenate dtseries within each (sub, ses) and write a global scan list."
    )
    ap.add_argument("--root", required=True,
                    help="Root directory to search recursively (BIDS-like tree).")
    ap.add_argument("--out", required=True,
                    help="Output directory for merged dtseries.")
    ap.add_argument("--pattern", default="*_desc-filtered_timeseries.dtseries.nii",
                    help="Glob pattern to find scans (default: %(default)s).")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print actions without executing.")
    ap.add_argument("--link-single", action="store_true",
                    help="For single-file sessions, create a symlink instead of copy.")
    ap.add_argument("--manifest", default=None,
                    help="Optional CSV manifest path (default: <out>/merge_manifest.csv).")
    ap.add_argument("--custom-order", nargs="*",
                    help="Optional ordering keywords (e.g., MID SST REST).")
    ap.add_argument("--scan-list-out", default=None,
                    help="Path to write 19_Scan_List_Concat.txt "
                         "(default: <pnet_inputs>/19_Scan_List_Concat.txt if <out> is inside pnet_inputs; "
                         "otherwise <out>/19_Scan_List_Concat.txt).")
    return ap.parse_args()

def extract_sub_ses(p: Path) -> Tuple[str, str]:
    sub = ses = None
    for part in p.parts:
        if not sub and BIDS_SUB_RE.fullmatch(part):
            sub = part
        if not ses and BIDS_SES_RE.fullmatch(part):
            ses = part
    return sub, ses

def sort_key(name: str, keywords: List[str]) -> Tuple[int, str]:
    lower = name.lower()
    rank = len(keywords)
    for i, kw in enumerate(keywords):
        if kw.lower() in lower:
            rank = i
            break
    return (rank, name)

def run(cmd: List[str], dry: bool):
    pretty = " ".join(f"'{c}'" if " " in c else c for c in cmd)
    print(">>", pretty)
    if not dry:
        subprocess.run(cmd, check=True)

def main():
    args = parse_args()
    wb = which("wb_command")
    root = Path(args.root).resolve()
    outdir = Path(args.out).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Decide default scan-list path
    if args.scan_list_out:
        scan_list_path = Path(args.scan_list_out).resolve()
    else:
        # If <out> is inside a directory named 'pnet_inputs', write list there; else write next to <out>
        parent = outdir.parent
        if parent.name == "pnet_inputs":
            scan_list_path = parent / "19_Scan_List_Concat.txt"
        else:
            scan_list_path = outdir / "19_Scan_List_Concat.txt"

    manifest_path = Path(args.manifest) if args.manifest else outdir / "merge_manifest.csv"

    print(f"[INFO] root={root}")
    print(f"[INFO] out ={outdir}")
    print(f"[INFO] pattern={args.pattern}")
    print(f"[INFO] dry-run={args.dry_run}  link-single={args.link_single}")
    print(f"[INFO] manifest={manifest_path}")
    print(f"[INFO] scan-list-out={scan_list_path}")
    if args.custom_order:
        print(f"[INFO] custom-order={args.custom_order}")

    # Discover files
    files = sorted(root.rglob(args.pattern))
    if not files:
        print(f"[WARN] No files found under {root} matching {args.pattern}")
        # still produce empty artifacts
        if not args.dry_run:
            manifest_path.parent.mkdir(parents=True, exist_ok=True)
            scan_list_path.parent.mkdir(parents=True, exist_ok=True)
            manifest_path.write_text("subject,session,n_inputs,output_path,inputs\n")
            scan_list_path.write_text("")
        print("[DONE] Nothing to merge.")
        return

    # Group by (sub, ses)
    groups: Dict[Tuple[str, str], List[Path]] = {}
    skipped = 0
    for f in files:
        sub, ses = extract_sub_ses(f)
        if not sub or not ses:
            skipped += 1
            continue
        groups.setdefault((sub, ses), []).append(f)
    if skipped:
        print(f"[WARN] Skipped {skipped} files without clear sub/ses in path.")

    merged_outputs: List[Path] = []
    rows = []

    # Process each group deterministically
    for (sub, ses) in sorted(groups.keys()):
        cands = groups[(sub, ses)]
        if args.custom_order:
            cands = sorted(cands, key=lambda p: sort_key(p.name, args.custom_order))
        else:
            cands = sorted(cands)

        n = len(cands)
        out = outdir / f"{sub}_{ses}_desc-merged_timeseries.dtseries.nii"

        if n == 0:
            print(f"[WARN] {sub} {ses}: 0 files — skipping.")
            rows.append([sub, ses, 0, "", ""])
            continue

        if n == 1:
            src = cands[0]
            print(f"[INFO] {sub} {ses}: 1 file → {'symlink' if args.link_single else 'copy'} → {out.name}")
            print(f"      - {src.name}")
            if not args.dry_run:
                if out.exists():
                    out.unlink()
                out.parent.mkdir(parents=True, exist_ok=True)
                if args.link_single:
                    rel = os.path.relpath(src, start=out.parent)
                    out.symlink_to(rel)
                else:
                    shutil.copy2(src, out)
            merged_outputs.append(out)
            rows.append([sub, ses, 1, str(out), str(src)])
            continue

        # n >= 2 → merge
        print(f"[INFO] {sub} {ses}: {n} files → merge → {out.name}")
        for p in cands:
            print(f"      - {p.name}")

        args_list = []
        for pth in cands:
            args_list.extend(["-cifti", str(pth)])

        cmd = [wb, "-cifti-merge", str(out)] + args_list
        run(cmd, args.dry_run)
        merged_outputs.append(out)
        rows.append([sub, ses, n, str(out), ";".join(str(p) for p in cands)])

    # Write manifest + scan list
    if not args.dry_run:
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        with open(manifest_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["subject", "session", "n_inputs", "output_path", "inputs"])
            w.writerows(rows)

        scan_list_path.parent.mkdir(parents=True, exist_ok=True)
        # Only include outputs that we created/would create
        with open(scan_list_path, "w") as f:
            for p in sorted({p.resolve() for p in merged_outputs}):
                f.write(str(p) + "\n")

    print(f"[DONE] Groups: {len(groups)}")
    print(f"[DONE] Manifest: {manifest_path}")
    if not args.dry_run:
        print(f"[DONE] Scan list: {scan_list_path}")

if __name__ == "__main__":
    main()
