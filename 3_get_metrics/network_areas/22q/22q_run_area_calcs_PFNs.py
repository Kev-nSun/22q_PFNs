from concurrent.futures import ProcessPoolExecutor, as_completed
import re
from pathlib import Path
import os, json
from area_calc_functions_22q import process_subject


surf_dir = Path("/cbica/projects/bbl_22q/data/derivatives_2025/22q_fsLR_32k_midthickness")
roi_dir  = Path("/cbica/projects/bbl_22q/analysis/allometry/inputs")
net_dir  = Path("/cbica/projects/bbl_22q/analysis/allometry/inputs/22q_PFN_loadings_normed")
deriv_dir  = Path("/cbica/projects/bbl_22q/analysis/allometry/outputs/22q_areas")
deriv_dir.mkdir(parents=True, exist_ok=True)

# Discover unique (sub, ses) from L-hemi files under sub-*/ses-*/anat/
pat = re.compile(
    r"^sub-(?P<sub>[^_]+)"
    r"_ses-(?P<ses>[^_]+)"
)

job_keys = []
for f in surf_dir.glob(
    "sub-*/ses-*/anat/sub-*_ses-*_hemi-L_space-fsLR_den-*_midthickness.surf.gii"
):
    m = pat.match(f.name)  # match filename only
    if not m:
        continue
    job_keys.append((m["sub"], m["ses"]))  # only keep relevant keys

job_keys = sorted(set(job_keys))

results, failures = [], {}
with ProcessPoolExecutor(max_workers=min(8, len(job_keys) or 1)) as ex:
    futs = {}  # map Future -> (sub, ses)
    for (sub, ses) in job_keys:
        sub_net_dir = net_dir / f"sub-{sub}"
        if not sub_net_dir.exists():  # Check PFNs exist by subject dir within net_dir
            failures[str((sub, ses))] = f"Missing net_dir: {sub_net_dir}"
            print(f"[SKIP] sub-{sub} ses-{ses}: PFN dir ({sub_net_dir}) not found")
            continue

        fut = ex.submit(
            process_subject,
            sub,
            ses,
            surf_dir,
            roi_dir,
            sub_net_dir,
            deriv_dir,
            atlas="PFN",
        )
        futs[fut] = (sub, ses)

    for fut in as_completed(futs):
        sub, ses = futs[fut]
        try:
            results.append(fut.result())
            print(f"[OK] sub-{sub} ses-{ses}")
        except Exception as e:
            failures[str((sub, ses))] = str(e)
            print(f"[FAIL] sub-{sub} ses-{ses}: {e}")

print(f"Ran {len(results)} subjects, skipped {len(failures)}.")

# Save summary CSV/JSON
(deriv_dir / "22q_summary_PFNs.json").write_text(json.dumps(results, indent=2))


import pandas as pd
REQUIRED = ("subject", "session", "TC_area", "network_areas")
rows = []
for r in results:
    missing = [k for k in REQUIRED if k not in r]
    if missing:
        raise KeyError(f"Result missing keys {missing}: {r}")
    row = {"subject": r["subject"], "session": r["session"], "TC_area": r["TC_area"]}
    row.update(r["network_areas"])
    rows.append(row)
df = pd.DataFrame(rows).set_index(["subject", "session"]).sort_index()
df.to_csv(deriv_dir/"22q_summary_PFNs.csv")