from __future__ import annotations
import subprocess, shutil, os
from pathlib import Path
import csv

def subject_dir(deriv: Path, sub: str, ses: str | None) -> Path:
    d = deriv / f"sub-{sub}"
    if ses:
        d = d / f"ses-{ses}"
    return d

def anat_dir(deriv: Path, sub: str, ses: str | None) -> Path:
    return subject_dir(deriv, sub, ses) / "anat"

def atlas_dir(deriv: Path, sub: str, ses: str | None, atlas: str) -> Path:
    return subject_dir(deriv, sub, ses) / atlas

def stats_dir(deriv: Path, sub: str, ses: str | None) -> Path:
    return subject_dir(deriv, sub, ses) / "stats"

def name_vertex_area_metric(sub, ses, den, hemi):
    ses_ent = f"_ses-{ses}" if ses else ""
    return f"sub-{sub}{ses_ent}_space-fsLR_den-{den}_hemi-{hemi}_vertex_area_metric.func.gii"

def name_vertex_area_dscalar(sub, ses, den):
    ses_ent = f"_ses-{ses}" if ses else ""
    return f"sub-{sub}{ses_ent}_space-fsLR_den-{den}_vertex_area_map.dscalar.nii"

def name_weighted_map(sub, ses, den, atlas, net):
    ses_ent = f"_ses-{ses}" if ses else ""
    return f"sub-{sub}{ses_ent}_space-fsLR_den-{den}_atlas-{atlas}_net-{net}_weightedarea_map.dscalar.nii"

def name_total_cortex_area_tsv(sub, ses, den):
    ses_ent = f"_ses-{ses}" if ses else ""
    return f"sub-{sub}{ses_ent}_space-fsLR_den-{den}_total_cortex_stat-SUM.tsv"

def name_network_areas_tsv(sub, ses, den, atlas):
    ses_ent = f"_ses-{ses}" if ses else ""
    return f"sub-{sub}{ses_ent}_space-fsLR_den-{den}_atlas-{atlas}_network_areas_stat-SUM.tsv"

def write_tsv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

def _resolve_wb(wb_command: str | None = None) -> str:
    if wb_command:
        return wb_command
    return shutil.which("wb_command") or os.environ.get("WB_COMMAND") or "/path/to/wb_command"

def run_wb(*args, capture: bool = False, wb_command: str | None = None) -> str:
    WB = _resolve_wb(wb_command)
    cmd = [WB, *map(str, args)]
    try:
        res = subprocess.run(
            cmd, check=True, text=True,
            stdout=subprocess.PIPE if capture else None,
            stderr=subprocess.PIPE
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"wb_command failed:\n{' '.join(cmd)}\nSTDERR:\n{e.stderr}") from e
    return res.stdout if capture else ""

def surface_vertex_areas(surf: Path, out_func: Path, wb_command: str | None = None):
    out_func.parent.mkdir(parents=True, exist_ok=True)
    run_wb("-surface-vertex-areas", surf, out_func, wb_command=wb_command)

def cifti_create_dense_scalar(out_dscalar: Path, left_metric: Path, right_metric: Path,
                              roi_left: Path, roi_right: Path, wb_command: str | None = None):
    run_wb("-cifti-create-dense-scalar", out_dscalar,
        "-left-metric", left_metric, "-roi-left", roi_left,
        "-right-metric", right_metric, "-roi-right", roi_right, 
        wb_command=wb_command)

def cifti_sum(in_dscalar: Path, wb_command: str | None = None) -> float:
    out = run_wb("-cifti-stats", in_dscalar, "-reduce", "SUM", capture=True, wb_command=wb_command)
    # wb_command prints a single number or a labeled line; be permissive:
    return float(out.strip().split()[-1])

def cifti_math(expr: str, out_dscalar: Path, wb_command: str | None = None, **vars_):
    args = ["-cifti-math", expr]
    for name, path in vars_.items():
        args += ["-var", name, Path(path)]
    args += [out_dscalar]
    run_wb(*args, wb_command=wb_command)

def process_subject(
    sub: str,
    surf_dir: Path,
    roi_dir: Path,
    net_dir: Path,
    deriv_dir: Path,
    roi_l: str | Path | None = None,
    roi_r: str | Path | None = None,
    wb_command: str | None = None,
    net_glob: str = "*.dscalar.nii",
    ses: str | None = "PNC1",
    acq: str = "refaced",
    density: str = "32k",
    atlas: str = "PNC_group",       # <- atlas label for filenames
):
    # Resolve BIDS-like directories
    ANAT = anat_dir(deriv_dir, sub, ses); ANAT.mkdir(parents=True, exist_ok=True)
    ATLS = atlas_dir(deriv_dir, sub, ses, atlas); ATLS.mkdir(parents=True, exist_ok=True)
    STATS = stats_dir(deriv_dir, sub, ses); STATS.mkdir(parents=True, exist_ok=True)

    # Input surface paths
    ses_ent = f"_ses-{ses}" if ses else ""
    surf_l = Path(surf_dir) / f"sub-{sub}{ses_ent}_acq-{acq}_hemi-L_space-fsLR_den-{density}_midthickness.surf.gii"
    surf_r = Path(surf_dir) / f"sub-{sub}{ses_ent}_acq-{acq}_hemi-R_space-fsLR_den-{density}_midthickness.surf.gii"

    # ROI paths
    roi_l = Path(roi_l) if roi_l is not None else Path(roi_dir) / f"S1200.L.atlasroi.{density}_fs_LR.shape.gii"
    roi_r = Path(roi_r) if roi_r is not None else Path(roi_dir) / f"S1200.R.atlasroi.{density}_fs_LR.shape.gii"

    # Check for inputs
    for p in (surf_l, surf_r, roi_l, roi_r):
        if not Path(p).exists():
            raise FileNotFoundError(f"Missing required file: {p}")

    # Outputs (BIDS-like names)
    area_l = ANAT / name_vertex_area_metric(sub, ses, density, "L")
    area_r = ANAT / name_vertex_area_metric(sub, ses, density, "R")
    area_cifti = ANAT / name_vertex_area_dscalar(sub, ses, density)

    # Compute vertex areas & combined dscalar
    if not area_l.exists(): surface_vertex_areas(surf_l, area_l, wb_command)  #get areas for L and R hemis
    if not area_r.exists(): surface_vertex_areas(surf_r, area_r, wb_command)
    if not area_cifti.exists():
        cifti_create_dense_scalar(area_cifti, area_l, area_r, roi_l, roi_r, wb_command) #combine L and R hemis into dscalar, excluding medial wall

    tc_area = cifti_sum(area_cifti, wb_command) #sum all per-vertex area values for total cortical (TC) area

    # Weighted networks
    weighted_sums = {}
    for net in sorted(Path(net_dir).glob(net_glob)):
        # derive network label for filename, e.g., PFN1_soft_parcel_normed -> PFN1
        # adjust regex to your actual filenames
        net_label = net.name.replace(".dscalar.nii", "")
        net_weighted_cifti = ATLS / name_weighted_map(sub, ses, density, atlas, net_label)
        if not net_weighted_cifti.exists():
            cifti_math("area * loading", net_weighted_cifti, wb_command, area=area_cifti, loading=net) # weight surface area values based on soft parcellation of each network
        weighted_sums[net_label] = cifti_sum(net_weighted_cifti, wb_command) # sum all per-vertex weighted area values for network area

    # Write subject-level stats TSVs
    write_tsv(STATS / name_total_cortex_area_tsv(sub, ses, density),
              [ {"subject": sub, "session": ses or "", "space": "fsLR", "den": density, "stat": "SUM", "TC_area": tc_area} ])

    rows = []
    for k, v in sorted(weighted_sums.items()):
        rows.append({"subject": sub, "session": ses or "", "space": "fsLR", "den": density,
                     "atlas": atlas, "network": k, "stat": "SUM", "area": v})
    write_tsv(STATS / name_network_areas_tsv(sub, ses, density, atlas), rows)

    return {"subject": sub, "TC_area": tc_area, "network_areas": weighted_sums}