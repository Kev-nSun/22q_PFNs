from concurrent.futures import ProcessPoolExecutor, as_completed
import re
from pathlib import Path
import os, json
from template_surf_area_calc_functions_PNC import process_subject


roi_dir = Path("/cbica/projects/bbl_22q/analysis/allometry/inputs")
net_dir = Path("/cbica/projects/bbl_22q/analysis/allometry/inputs/PNC_PFN_loadings_normed")
deriv_dir = Path("/cbica/projects/bbl_22q/analysis/allometry/outputs/PNC_areas")
deriv_dir.mkdir(parents=True, exist_ok=True)

ses = "PNC1"
density = "32k"
atlas = "PFN"

# Discover subjects from available PFN-loading directories.
job_keys = []
for d in sorted(net_dir.glob("sub-*")):
    if not d.is_dir():
        continue

    sub = d.name.replace("sub-", "")

    # Require at least one network loading file.
    if not any(d.glob("*.dscalar.nii")):
        print(f"[SKIP] sub-{sub}: no PFN loading files found in {d}")
        continue

    job_keys.append((sub, ses, density))

results, failures = [], {}

max_workers = min(8, len(job_keys)) if job_keys else 1

with ProcessPoolExecutor(max_workers=max_workers) as ex:
    futs = {}

    for sub, ses_i, den in job_keys:
        sub_net_dir = net_dir / f"sub-{sub}"

        fut = ex.submit(
            process_subject,
            sub,
            roi_dir,
            sub_net_dir,
            deriv_dir,
            ses=ses_i,
            density=den,
            atlas=atlas,
        )

        futs[fut] = (sub, ses_i, den)

    for fut in as_completed(futs):
        sub, ses_i, den = futs[fut]

        try:
            results.append(fut.result())
            print(f"[OK] sub-{sub} ses-{ses_i} den-{den}")
        except Exception as e:
            failures[str((sub, ses_i, den))] = str(e)
            print(f"[FAIL] sub-{sub} ses-{ses_i} den-{den}: {e}")

print(f"Ran {len(results)} subjects, failed/skipped {len(failures)}.")

summary = {
    "results": results,
    "failures": failures,
}

(deriv_dir / "PNC_template_summary_PFNs.json").write_text(json.dumps(summary, indent=2))

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
        df = pd.DataFrame(rows).set_index("subject").sort_index()
        df.to_csv(deriv_dir / "PNC_template_summary_PFNs.csv")

except Exception as e:
    print(f"[WARN] Could not write CSV summary: {e}")