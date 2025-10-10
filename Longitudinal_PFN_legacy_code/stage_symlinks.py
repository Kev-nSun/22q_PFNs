#python stage_symlinks.py \
#  --root /cbica/projects/PFN_ABCD/abcd-hcp-pipeline_0.1.4_timeseries \
#  --stage-out /cbica/projects/PFN_ABCD/pnet_inputs/pnet_stage_symlinks \
#  --scans-out /cbica/projects/PFN_ABCD/pnet_inputs/pnet_stage_symlinks/19_Scan_List.txt

#!/usr/bin/env python3
import argparse
import re
import sys
from pathlib import Path

BIDS_NAME_RE = re.compile(r"(sub-[A-Za-z0-9]+)_(ses-[A-Za-z0-9]+)")

def parse_sub_ses(p: Path):
    """Try filename first; fall back to path components."""
    m = BIDS_NAME_RE.search(p.name)
    if m:
        return m.group(1), m.group(2)
    # fallback: search any path part
    parts = "/".join(p.parts)
    m2 = BIDS_NAME_RE.search(parts)
    if m2:
        return m2.group(1), m2.group(2)
    return None, None

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Stage *_desc-filtered_timeseries.dtseries.nii into a symlink tree grouped by sub+ses.\n"
            "Each group is a folder named 'sub-XXXX+ses-YYYY' containing symlinks to its scans.\n"
            "Optionally writes a file_scans.txt listing the symlink paths (one per line)."
        )
    )
    ap.add_argument("--root", required=True,
                    help="Root directory to search recursively (your timeseries tree)")
    ap.add_argument("--stage-out", required=True,
                    help="Directory to create the staged symlink tree (will be created if missing)")
    ap.add_argument("--pattern", default="*_desc-filtered_timeseries.dtseries.nii",
                    help="Glob pattern to find scans (default: %(default)s)")
    ap.add_argument("--scans-out", default=None,
                    help="Optional path to write a file_scans.txt listing the **symlink** paths")
    ap.add_argument("--overwrite", action="store_true",
                    help="Replace existing symlinks if they already exist")
    ap.add_argument("--dry-run", action="store_true",
                    help="Show what would be done without creating anything")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    stage = Path(args.stage_out).resolve()
    if not root.exists():
        sys.exit(f"[ERROR] --root does not exist: {root}")

    # Collect candidate files
    files = sorted(root.rglob(args.pattern))
    if not files:
        sys.exit(f"[ERROR] Found 0 files under {root} matching pattern {args.pattern}")

    # Build mapping: key -> list[files]
    groups = {}  # key: "sub+ses" -> list[Path]
    skipped = []
    for f in files:
        sub, ses = parse_sub_ses(f)
        if not sub or not ses:
            skipped.append(f)
            continue
        key = f"{sub}+{ses}"
        groups.setdefault(key, []).append(f)

    if skipped:
        print(f"[WARN] Skipped {len(skipped)} files without sub/ses in name or path.")
        for s in skipped[:5]:
            print(f"       - {s}")
        if len(skipped) > 5:
            print("       ...")

    # Prepare stage dir
    if args.dry_run:
        print(f"[DRY-RUN] Would create stage root: {stage}")
    else:
        stage.mkdir(parents=True, exist_ok=True)

    # Create symlinks, deterministic order, numbered names
    symlink_paths = []
    for key in sorted(groups.keys()):
        dest_dir = stage / key
        if args.dry_run:
            print(f"[DRY-RUN] Would create dir: {dest_dir}")
        else:
            dest_dir.mkdir(exist_ok=True)

        scans = sorted(groups[key])  # deterministic
        for index, source in enumerate(scans, start=1):
            # Keep original basename for traceability, but prefix with an index
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
                        # already there; keep it
                        symlink_paths.append(dest)
                        continue
                # Create relative symlink for portability
                rel = Path(Path.cwd()).joinpath(source).resolve()
                dest.symlink_to(rel)
            symlink_paths.append(dest)

    # Optionally write file_scans.txt listing the symlinked paths
    if args.scans_out:
        scans_out = Path(args.scans_out).resolve()
        if args.dry_run:
            print(f"[DRY-RUN] Would write {len(symlink_paths)} paths to {scans_out}")
        else:
            scans_out.parent.mkdir(parents=True, exist_ok=True)
            scans_out.write_text("\n".join(str(p) for p in symlink_paths))
            print(f"[OK] Wrote scan list: {scans_out} ({len(symlink_paths)} entries)")

    # Summary
    total = sum(len(v) for v in groups.values())
    print(f"[DONE] Groups: {len(groups)}  |  Symlinks: {total}  |  Stage root: {stage}")

    print("\nNext steps:")
    print("  • Set pNet config to use this stage for automatic concatenation>")
    print("      necessary_settings.file_scans = '<path to your scans_out>'  (recommended)")
    print("      OR point to all symlink files if you don’t write scans_out.")
    print("  • Do NOT set file_subject_ID or file_subject_folder.")
    print("  • Ensure Combine_Scan = true.")

if __name__ == "__main__":
    main()