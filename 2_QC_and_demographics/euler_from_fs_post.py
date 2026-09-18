#This script wrties a summary Euler csv of 22q subjects based on tsv QC files produced by freesurfer-post
#Run: python euler_from_fs_post.py /cbica/projects/bbl_22q/data/derivatives_2025/22q_freesurfer-post_QC/22q_freesurfer_qc_euler_summary.csv
#OR Run: python euler_from_fs_post.py /cbica/projects/bbl_22q/data/derivatives_2025/PNC_freesurfer-post_QC/PNC_freesurfer_qc_euler_summary.csv

import csv
import glob
import sys

#PATTERN = "/cbica/projects/bbl_22q/data/derivatives_2025/22q_freesurfer-post_QC/freesurfer-post/sub-*/sub-*_ses-*_desc-FreeSurfer_qc.tsv"
PATTERN = "/cbica/projects/bbl_22q/data/derivatives_2025/PNC_freesurfer-post_QC/sub-*_ses-PNC1_desc-FreeSurfer_qc.tsv"

def find_col_indices(header):
    """
    Returns indices for participant_id, session_id, lh_euler, rh_euler.
    If participant_id/session_id appear multiple times, uses the first occurrence.
    """
    def first_idx(name):
        for i, h in enumerate(header):
            if h == name:
                return i
        return None

    idx = {
        "participant_id": first_idx("participant_id"),
        "session_id": first_idx("session_id"),
        "lh_euler": first_idx("lh_euler"),
        "rh_euler": first_idx("rh_euler"),
    }
    missing = [k for k, v in idx.items() if v is None]
    if missing:
        raise ValueError(f"Missing required columns in header: {missing}")
    return idx

def coerce_float(x):
    s = str(x).strip()
    if s == "" or s.lower() in {"na", "nan", "none"}:
        return None
    return float(s)

def read_one_row(tsv_path):
    with open(tsv_path, "r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if header is None:
            raise ValueError("Empty file")

        idx = find_col_indices(header)

        row = next(reader, None)
        if row is None:
            raise ValueError("No data rows found")

        sub = row[idx["participant_id"]].strip()
        ses = row[idx["session_id"]].strip()
        lh = coerce_float(row[idx["lh_euler"]])
        rh = coerce_float(row[idx["rh_euler"]])
        euler = (lh + rh) if (lh is not None and rh is not None) else None

        return sub, ses, lh, rh, euler

def main(out_csv):
    tsvs = sorted(glob.glob(PATTERN))
    if not tsvs:
        print(f"No files matched pattern:\n  {PATTERN}", file=sys.stderr)
        return 2

    rows = []
    for tsv in tsvs:
        try:
            sub, ses, lh, rh, euler = read_one_row(tsv)
            rows.append({
                "subject": sub,
                "session": ses,
                "lh_euler": lh,
                "rh_euler": rh,
                "euler": euler
            })
        except Exception as e:
            print(f"ERROR: {tsv}: {e}", file=sys.stderr)

    # Write output CSV
    with open(out_csv, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=["subject", "session", "lh_euler", "rh_euler", "euler"])
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} rows -> {out_csv}")
    return 0

if __name__ == "__main__":
    #out_csv = sys.argv[1] if len(sys.argv) > 1 else "/cbica/projects/bbl_22q/data/derivatives_2025/PNC_freesurfer-post_QC/PNC_freesurfer_qc_euler_summary.csv"
    out_csv = sys.argv[1] if len(sys.argv) > 1 else "/cbica/projects/bbl_22q/data/derivatives_2025/22q_freesurfer-post_QC/22q_freesurfer_qc_euler_summary.csv"
    raise SystemExit(main(out_csv))
