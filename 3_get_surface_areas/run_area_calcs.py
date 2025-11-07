from concurrent.futures import ProcessPoolExecutor, as_completed
import re
from pathlib import Path
import os, json
from area_calc_functions import process_subject


surf_dir = Path("/cbica/projects/bbl_22q/data/PNC_fsLR_32k_midthickness")
roi_dir  = Path("/cbica/projects/bbl_22q/analysis/allometry/inputs")
net_dir  = Path("/cbica/projects/bbl_22q/analysis/allometry/inputs/PNC_group_atlas_normed")
deriv_dir  = Path("/cbica/projects/bbl_22q/analysis/allometry/outputs")
deriv_dir.mkdir(parents=True, exist_ok=True)

# Discover unique (sub, ses, acq, den) from L-hemi files
pat = re.compile(
    r"^sub-(?P<sub>[^_]+)"
    r"(?:_ses-(?P<ses>[^_]+))?"
    r"(?:_acq-(?P<acq>[^_]+))?"
    r"_hemi-L_space-fsLR_den-(?P<den>[^_]+)_midthickness\.surf\.gii$"
)
job_keys = []
for f in surf_dir.glob("sub-*_hemi-L_space-fsLR_den-*_midthickness.surf.gii"):
    m = pat.match(f.name)
    if not m:
        continue
    job_keys.append( (m["sub"], m["ses"] or "PNC1", m["acq"] or "refaced", m["den"] or "32k") )
job_keys = sorted(set(job_keys))

results, failures = [], {}
with ProcessPoolExecutor(max_workers=min(8, len(job_keys))) as ex:
    futs = {
        ex.submit(process_subject, sub, surf_dir, roi_dir, net_dir, deriv_dir,
                  ses=ses, acq=acq, density=den, atlas="PNC_group"): (sub, ses, acq, den)
        for (sub, ses, acq, den) in job_keys
    }
    for fut in as_completed(futs):
        key = futs[fut]
        try:
            results.append(fut.result())
            print(f"[OK] sub-{key[0]} ses-{key[1]} den-{key[3]}")
        except Exception as e:
            failures[str(key)] = str(e)
            print(f"[FAIL] sub-{key}: {e}")

# Save summary CSV/JSON
(deriv_dir / "summary.json").write_text(json.dumps(results, indent=2))

try:
    import pandas as pd
    rows = []
    for r in results:
        row = {"subject": r["subject"], "TC_area": r["TC_area"]}
        row.update(r["network_areas"])
        rows.append(row)
    df = pd.DataFrame(rows).set_index("subject").sort_index()
    df.to_csv(deriv_dir/"summary.csv")
except Exception:
    pass