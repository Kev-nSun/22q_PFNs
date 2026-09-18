from __future__ import annotations

import subprocess
import shutil
import os
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


def name_midthickness_surface(sub, ses, den, hemi):
    ses_ent = f"_ses-{ses}" if ses else ""
    return (
        f"sub-{sub}{ses_ent}_space-fsLR_den-{den}_hemi-{hemi}"
        "_desc-template_midthickness.surf.gii"
    )


def name_vertex_area_metric(sub, ses, den, hemi):
    ses_ent = f"_ses-{ses}" if ses else ""
    return (
        f"sub-{sub}{ses_ent}_space-fsLR_den-{den}_hemi-{hemi}"
        "_desc-template_vertexarea.shape.gii"
    )


def name_vertex_area_dscalar(sub, ses, den):
    ses_ent = f"_ses-{ses}" if ses else ""
    return (
        f"sub-{sub}{ses_ent}_space-fsLR_den-{den}"
        "_desc-template_vertexarea.dscalar.nii"
    )


def name_weighted_map(sub, ses, den, atlas, net):
    ses_ent = f"_ses-{ses}" if ses else ""
    return (
        f"sub-{sub}{ses_ent}_space-fsLR_den-{den}_atlas-{atlas}"
        f"_net-{net}_desc-template_weightedarea_map.dscalar.nii"
    )


def name_total_cortex_area_tsv(sub, ses, den):
    ses_ent = f"_ses-{ses}" if ses else ""
    return (
        f"sub-{sub}{ses_ent}_space-fsLR_den-{den}"
        "_desc-template_total_cortex_stat-SUM.tsv"
    )


def name_network_areas_tsv(sub, ses, den, atlas):
    ses_ent = f"_ses-{ses}" if ses else ""
    return (
        f"sub-{sub}{ses_ent}_space-fsLR_den-{den}_atlas-{atlas}"
        "_desc-template_network_areas_stat-SUM.tsv"
    )


def write_tsv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)

    if not rows:
        return

    fieldnames = list(rows[0].keys())

    with path.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=fieldnames,
        )
        w.writeheader()
        w.writerows(rows)


def _resolve_wb(wb_command: str | None = None) -> str:
    if wb_command:
        return wb_command

    return (
        shutil.which("wb_command")
        or os.environ.get("WB_COMMAND")
        or "/path/to/wb_command"
    )


def run_wb(
    *args,
    capture: bool = False,
    wb_command: str | None = None,
) -> str:

    WB = _resolve_wb(wb_command)

    cmd = [WB, *map(str, args)]

    try:
        res = subprocess.run(
            cmd,
            check=True,
            text=True,
            stdout=subprocess.PIPE if capture else None,
            stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"wb_command failed:\n"
            f"{' '.join(cmd)}\n"
            f"STDERR:\n{e.stderr}"
        ) from e

    return res.stdout if capture else ""


def surface_vertex_areas(
    surface: Path,
    out_metric: Path,
    wb_command: str | None = None,
):
    """
    Compute geometric vertex areas directly from a surface mesh.

    Workbench assigns each vertex one third of the area of every
    triangle incident on that vertex.

    Output units: mm^2.
    """
    out_metric.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    run_wb(
        "-surface-vertex-areas",
        surface,
        out_metric,
        wb_command=wb_command,
    )


def cifti_create_dense_scalar(
    out_dscalar: Path,
    left_metric: Path,
    right_metric: Path,
    roi_left: Path,
    roi_right: Path,
    wb_command: str | None = None,
):

    out_dscalar.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    run_wb(
        "-cifti-create-dense-scalar",
        out_dscalar,
        "-left-metric",
        left_metric,
        "-roi-left",
        roi_left,
        "-right-metric",
        right_metric,
        "-roi-right",
        roi_right,
        wb_command=wb_command,
    )


def cifti_sum(
    in_dscalar: Path,
    wb_command: str | None = None,
) -> float:

    out = run_wb(
        "-cifti-stats",
        in_dscalar,
        "-reduce",
        "SUM",
        capture=True,
        wb_command=wb_command,
    )

    return float(out.strip().split()[-1])


def cifti_math(
    expr: str,
    out_dscalar: Path,
    wb_command: str | None = None,
    **vars_,
):

    out_dscalar.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    args = [
        "-cifti-math",
        expr,
        out_dscalar,
    ]

    for name, path in vars_.items():
        args += [
            "-var",
            name,
            Path(path),
        ]

    run_wb(
        *args,
        wb_command=wb_command,
    )


def _link_or_copy(
    src: Path,
    dst: Path,
    use_symlink: bool = True,
):

    dst.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    if dst.exists() or dst.is_symlink():
        return

    if use_symlink:
        try:
            dst.symlink_to(src)
            return
        except OSError:
            pass

    shutil.copy2(src, dst)


def process_subject(
    sub: str,
    roi_dir: Path,
    net_dir: Path,
    deriv_dir: Path,
    roi_l: str | Path | None = None,
    roi_r: str | Path | None = None,
    wb_command: str | None = None,
    net_glob: str = "*.dscalar.nii",
    ses: str | None = "PNC1",
    density: str = "32k",
    atlas: str = "PNC_group",
    use_symlink: bool = True,
):
    """
    Compute template-fsLR surface-area-weighted network areas.

    IMPORTANT:
    Vertex areas are calculated directly from the geometric
    TemplateFlow fsLR midthickness surface using:

        wb_command -surface-vertex-areas

    Therefore these are the literal geometric vertex areas of the
    distributed fsLR template midthickness mesh.

    They are NOT TemplateFlow desc-vaavg group-average vertex areas.

    Results therefore represent network loading integrated over the
    geometric fsLR template surface.
    """

    ANAT = anat_dir(
        deriv_dir,
        sub,
        ses,
    )

    ATLS = atlas_dir(
        deriv_dir,
        sub,
        ses,
        atlas,
    )

    STATS = stats_dir(
        deriv_dir,
        sub,
        ses,
    )

    ANAT.mkdir(
        parents=True,
        exist_ok=True,
    )

    ATLS.mkdir(
        parents=True,
        exist_ok=True,
    )

    STATS.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ------------------------------------------------------------
    # Cortical ROI
    # ------------------------------------------------------------

    roi_l = (
        Path(roi_l)
        if roi_l is not None
        else Path(roi_dir)
        / f"S1200.L.atlasroi.{density}_fs_LR.shape.gii"
    )

    roi_r = (
        Path(roi_r)
        if roi_r is not None
        else Path(roi_dir)
        / f"S1200.R.atlasroi.{density}_fs_LR.shape.gii"
    )

    # ------------------------------------------------------------
    # Get the actual TemplateFlow fsLR midthickness surfaces
    # ------------------------------------------------------------

    template_surf_dir = Path(
        "/cbica/projects/bbl_22q/templateflow_home/tpl-fsLR"
    )

    tf_surf_l = (
        template_surf_dir
        / "tpl-fsLR_den-32k_hemi-L_midthickness.surf.gii"
    )

    tf_surf_r = (
        template_surf_dir
        / "tpl-fsLR_den-32k_hemi-R_midthickness.surf.gii"
    )

    for p in (
        tf_surf_l,
        tf_surf_r,
        roi_l,
        roi_r,
    ):
        if not Path(p).exists():
            raise FileNotFoundError(
                f"Missing required file: {p}"
            )

    # ------------------------------------------------------------
    # Put template surfaces in derivatives directory
    # ------------------------------------------------------------

    surf_l = ANAT / name_midthickness_surface(
        sub,
        ses,
        density,
        "L",
    )

    surf_r = ANAT / name_midthickness_surface(
        sub,
        ses,
        density,
        "R",
    )

    _link_or_copy(
        tf_surf_l,
        surf_l,
        use_symlink=use_symlink,
    )

    _link_or_copy(
        tf_surf_r,
        surf_r,
        use_symlink=use_symlink,
    )

    # ------------------------------------------------------------
    # Calculate ACTUAL geometric template vertex areas
    # ------------------------------------------------------------

    area_l = ANAT / name_vertex_area_metric(
        sub,
        ses,
        density,
        "L",
    )

    area_r = ANAT / name_vertex_area_metric(
        sub,
        ses,
        density,
        "R",
    )

    if not area_l.exists():
        surface_vertex_areas(
            surf_l,
            area_l,
            wb_command=wb_command,
        )

    if not area_r.exists():
        surface_vertex_areas(
            surf_r,
            area_r,
            wb_command=wb_command,
        )

    # ------------------------------------------------------------
    # Combine hemisphere metrics into cortical CIFTI
    # ------------------------------------------------------------

    area_cifti = ANAT / name_vertex_area_dscalar(
        sub,
        ses,
        density,
    )

    if not area_cifti.exists():
        cifti_create_dense_scalar(
            area_cifti,
            area_l,
            area_r,
            roi_l,
            roi_r,
            wb_command=wb_command,
        )

    # Total geometric area of the template cortex,
    # restricted by the cortical ROI masks.
    tc_area = cifti_sum(
        area_cifti,
        wb_command,
    )

    # ------------------------------------------------------------
    # Surface-area-weight each network map
    # ------------------------------------------------------------

    weighted_sums = {}

    for net in sorted(
        Path(net_dir).glob(net_glob)
    ):

        net_label = net.name.replace(
            ".dscalar.nii",
            "",
        )

        net_weighted_cifti = (
            ATLS
            / name_weighted_map(
                sub,
                ses,
                density,
                atlas,
                net_label,
            )
        )

        if not net_weighted_cifti.exists():
            cifti_math(
                "area * loading",
                net_weighted_cifti,
                wb_command,
                area=area_cifti,
                loading=net,
            )

        weighted_sums[net_label] = cifti_sum(
            net_weighted_cifti,
            wb_command,
        )

    # ------------------------------------------------------------
    # Total cortex area TSV
    # ------------------------------------------------------------

    write_tsv(
        STATS
        / name_total_cortex_area_tsv(
            sub,
            ses,
            density,
        ),
        [
            {
                "subject": sub,
                "session": ses or "",
                "space": "fsLR",
                "den": density,
                "surface": "template_midthickness",
                "area_source": (
                    "TemplateFlow_midthickness_geometry"
                ),
                "stat": "SUM",
                "TC_area": tc_area,
            }
        ],
    )

    # ------------------------------------------------------------
    # Network area TSV
    # ------------------------------------------------------------

    rows = []

    for network, area in sorted(
        weighted_sums.items()
    ):
        rows.append(
            {
                "subject": sub,
                "session": ses or "",
                "space": "fsLR",
                "den": density,
                "surface": "template_midthickness",
                "area_source": (
                    "TemplateFlow_midthickness_geometry"
                ),
                "atlas": atlas,
                "network": network,
                "stat": "SUM",
                "area": area,
            }
        )

    write_tsv(
        STATS
        / name_network_areas_tsv(
            sub,
            ses,
            density,
            atlas,
        ),
        rows,
    )

    return {
        "subject": sub,
        "session": ses or "",
        "surface": "template_midthickness",
        "area_source": (
            "TemplateFlow_midthickness_geometry"
        ),
        "TC_area": tc_area,
        "network_areas": weighted_sums,
    }