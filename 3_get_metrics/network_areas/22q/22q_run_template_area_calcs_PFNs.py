from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import json
import re

from template_surf_area_calc_functions_22q import process_subject


roi_dir = Path("/cbica/projects/bbl_22q/analysis/allometry/inputs")
net_dir = Path("/cbica/projects/bbl_22q/analysis/allometry/inputs/22q_PFN_loadings_normed")
deriv_dir = Path("/cbica/projects/bbl_22q/analysis/allometry/outputs/22q_areas")
deriv_dir.mkdir(parents=True, exist_ok=True)

density = "32k"
atlas = "PFN"

ses_pat = re.compile(r"_ses-([^_]+)")

job_keys = []

# Discover jobs from available PFN-loading folders, not midthickness surfaces.
for sub_dir in sorted(net_dir.glob("sub-*")):
    if not sub_dir.is_dir():
        continue

    sub = sub_dir.name.replace("sub-", "")

    # Case 1: subject/session nested PFN dirs: sub-*/ses-*/*.dscalar.nii
    ses_dirs = sorted([d for d in sub_dir.glob("ses-*") if d.is_dir()])
    if ses_dirs:
        for ses_dir in ses_dirs:
            ses = ses_dir.name.replace("ses-", "")
            if any(ses_dir.glob("*.dscalar.nii")):
                job_keys.append((sub, ses, ses_dir))
            else:
                print(f"[SKIP] sub-{sub} ses-{ses}: no PFN loading files in {ses_dir}")
        continue

    # Case 2: flat subject PFN dir: sub-*/*.dscalar.nii
    files = sorted(sub_dir.glob("*.dscalar.nii"))
    if not files:
        print(f"[SKIP] sub-{sub}: no PFN loading files in {sub_dir}")
        continue

    # Infer sessions from filenames if present.
    sessions = sorted(
        {
            m.group(1)
            for f in files
            if (m := ses_pat.search(f.name))
        }
    )

    if sessions:
        for ses in sessions:
            job_keys.append((sub, ses, sub_dir))
    else:
        # Fallback if loading files do not encode session.
        # Change this default if your 22q data use a different session label.
        job_keys.append((sub, "1", sub_dir))

job_keys = sorted(set(job_keys))

results, failures = [], {}

with ProcessPoolExecutor(max_workers=min(8, len(job_keys) or 1)) as ex:
    futs = {}

    for sub, ses, sub_net_dir in job_keys:
        fut = ex.submit(
            process_subject,
            sub,
            ses,
            roi_dir,
            sub_net_dir,
            deriv_dir,
            density=density,
            atlas=atlas,
        )
        futs[fut] = (sub, ses, sub_net_dir)

    for fut in as_completed(futs):
        sub, ses, sub_net_dir = futs[fut]

        try:
            results.append(fut.result())
            print(f"[OK] sub-{sub} ses-{ses}")
        except Exception as e:
            failures[str((sub, ses))] = str(e)
            print(f"[FAIL] sub-{sub} ses-{ses}: {e}")

print(f"Ran {len(results)} subjects/sessions, failed/skipped {len(failures)}.")

summary = {
    "results": results,
    "failures": failures,
}

(deriv_dir / "22q_template_summary_PFNs.json").write_text(json.dumps(summary, indent=2))

try:
    import pandas as pd

    rows = []
    for r in results:
        row = {
            "subject": r["subject"],
            "session": r.get("session", ""),
            "surface": r.get("surface", ""),
            "area_source": r.get("area_source", ""),
            "TC_area": r["TC_area"],
        }
        row.update(r["network_areas"])
        rows.append(row)

    if rows:
        df = pd.DataFrame(rows).set_index(["subject", "session"]).sort_index()
        df.to_csv(deriv_dir / "22q_template_summary_PFNs.csv")

except Exception as e:
    print(f"[WARN] Could not write CSV summary: {e}")