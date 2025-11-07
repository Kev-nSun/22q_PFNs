#python stage_symlinks.py \
#  --root /cbica/projects/bbl_22q/data/22q_xcpd_extracted \
#  --stage-out /cbica/projects/bbl_22q/analysis/pfn/pnet_inputs/22q_pnet_stage_symlinks \
#  --scans-out /cbica/projects/bbl_22q/analysis/pfn/pnet_inputs/22q_pnet_stage_symlinks/Scan_List.txt \
#   --group subject

import argparse
import re
import sys
import os
from pathlib import Path

# Match sub and optional ses anywhere in the filename/path
SUB_RE = re.compile(r"(sub-[A-Za-z0-9]+)")
SES_RE = re.compile(r"(ses-[A-Za-z0-9]+)")

def parse_sub_ses(p: Path):
    """Return (sub, ses or None). Search filename first; then full path."""
    name = p.name
    sub = SUB_RE.search(name) or SUB_RE.search("/".join(p.parts))
    ses = SES_RE.search(name) or SES_RE.search("/".join(p.parts))
    sub = sub.group(1) if sub else None
    ses = ses.group(1) if ses else None
    return sub, ses

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Stage *.dtseries.nii into a symlink tree grouped by subject (or subject+session).\n"
            "Each group is a folder containing symlinks to its scans.\n"
            "Optionally writes a file_scans.txt listing the **symlink** paths (one per line)."
        )
    )
    ap.add_argument("--root", required=True,
                    help="Root directory to search recursively (your timeseries tree)")
    ap.add_argument("--stage-out", required=True,
                    help="Directory to create the staged symlink tree (will be created if missing)")
    ap.add_argument("--pattern", default="*.dtseries.nii",
                    help="Glob pattern to find scans (default: %(default)s)")
    ap.add_argument("--scans-out", default=None,
                    help="Optional path to write a file_scans.txt listing the **symlink** paths")
    ap.add_argument("--overwrite", action="store_true",
                    help="Replace existing symlinks if they already exist")
    ap.add_argument("--dry-run", action="store_true",
                    help="Show what would be done without creating anything")
    ap.add_argument("--group", choices=["subject", "subject+session"], default="subject+session",
                    help="Grouping key for concatenation (default: subject+session)")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    stage = Path(args.stage_out).resolve()
    if not root.exists():
        sys.exit(f"[ERROR] --root does not exist: {root}")

    files = sorted(root.rglob(args.pattern))
    if not files:
        sys.exit(f"[ERROR] Found 0 files under {root} matching pattern {args.pattern}")

    groups = {}  # key -> list[Path]
    skipped = []
    for f in files:
        sub, ses = parse_sub_ses(f)
        if not sub:
            skipped.append(f)
            continue
        if args.group == "subject":
            key = f"{sub}"
        else:  # subject+session
            if not ses:
                skipped.append(f)
                continue
            key = f"{sub}+{ses}"
        groups.setdefault(key, []).append(f)

    if skipped:
        print(f"[WARN] Skipped {len(skipped)} files without required IDs.")
        for s in skipped[:5]:
            print(f"       - {s}")
        if len(skipped) > 5:
            print("       ...")

    if args.dry_run:
        print(f"[DRY-RUN] Would create stage root: {stage}")
    else:
        stage.mkdir(parents=True, exist_ok=True)

    symlink_paths = []
    for key in sorted(groups.keys()):
        dest_dir = stage / key
        if args.dry_run:
            print(f"[DRY-RUN] Would create dir: {dest_dir}")
        else:
            dest_dir.mkdir(exist_ok=True)

        scans = sorted(groups[key])
        for index, source in enumerate(scans, start=1):
            out_name = f"{index:03d}__{source.name}"
            dest = dest_dir / out_name

            if args.dry_run:
                action = "overwrite" if (dest.exists() and args.overwrite) else "create"
                print(f"[DRY-RUN] Would {action} symlink: {dest} -> {source}")
            else:
                if dest.exists():
                    if args.overwrite:
                        dest.unlink()
                    else:
                        symlink_paths.append(dest)
                        continue
                # Create a true relative symlink (portable if stage folder moves with tree)
                rel = os.path.relpath(source.resolve(), start=dest.parent.resolve())
                dest.symlink_to(rel)
            symlink_paths.append(dest)

    if args.scans_out:
        scans_out = Path(args.scans_out).resolve()
        if args.dry_run:
            print(f"[DRY-RUN] Would write {len(symlink_paths)} paths to {scans_out}")
        else:
            scans_out.parent.mkdir(parents=True, exist_ok=True)
            scans_out.write_text("\n".join(str(p) for p in symlink_paths))
            print(f"[OK] Wrote scan list: {scans_out} ({len(symlink_paths)} entries)")

    total = sum(len(v) for v in groups.values())
    print(f"[DONE] Groups: {len(groups)}  |  Symlinks: {total}  |  Stage root: {stage}")

    print("\nNext steps:")
    print("  • Set pNet config to use this stage for automatic concatenation>")
    print("  • Do NOT set file_subject_ID or file_subject_folder.")
    print("  • Ensure Combine_Scan = "True".")

if __name__ == "__main__":
    main()