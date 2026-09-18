#!/usr/bin/env python3

import subprocess
from pathlib import Path
import re

# -------- user settings --------
wb_command = "wb_command"

base_dir = Path(
    "C:/Users/kevin/OneDrive/Documents/NGG_PhD/"
    "Alexander-Bloch/22q_Project/Analyses/Results/"
    "TCA_PFN_loadings_gam_predictions/091426_VisualizeFolder"
)

scene_dir = base_dir / "100_frame_Scene_dir"
scene_dir.mkdir(parents=True, exist_ok=True)

constant_dlabel = Path("C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/PNC_data/PNC_group_atlas/PNC_group_hard_parcellation_mwall.dlabel.nii")
# --------------------------------


def extract_frame_number(path: Path) -> int:

    pattern = r"Frame_(\d+)_mwall\.dlabel\.nii$"

    m = re.match(pattern, path.name)

    if not m:
        raise ValueError(
            f"Could not extract frame number from filename: {path.name}"
        )

    return int(m.group(1))


def check_frames(found_frames):

    expected = list(range(1, 101))

    if found_frames != expected:
        raise RuntimeError(
            f"Frame mismatch! Expected 1..100, found: {found_frames}"
        )

    return expected


# ============================================================
# DYNAMIC PFN FILE
# ============================================================

dynamic_output = (
    scene_dir /
    "PFN_loading_normed_movie_frames_1-99q_mwall.dlabel.nii"
)

input_dlabels = list(
    base_dir.glob("Frame_*/Frame_*_mwall.dlabel.nii")
)

if not input_dlabels:
    raise RuntimeError(
        f"No dlabel files found in {base_dir}"
    )

input_dlabels = sorted(
    input_dlabels,
    key=extract_frame_number
)

found_frames = [
    extract_frame_number(p)
    for p in input_dlabels
]

expected_frames = check_frames(found_frames)


dynamic_cmd = [
    wb_command,
    "-cifti-merge",
    str(dynamic_output)
]

for f in input_dlabels:
    dynamic_cmd.extend([
        "-cifti",
        str(f),
        "-index",
        "1"
    ])


print("\n========================================")
print("Creating dynamic 100-frame dlabel")
print("========================================")
print(f"Output: {dynamic_output}")

subprocess.run(dynamic_cmd, check=True)


# ============================================================
# CONSTANT OVERLAY FILE
# ============================================================

constant_output = (
    scene_dir /
    "group_atlas_100frames.dlabel.nii"
)

if not constant_dlabel.exists():
    raise RuntimeError(
        f"Constant dlabel does not exist:\n{constant_dlabel}"
    )


constant_cmd = [
    wb_command,
    "-cifti-merge",
    str(constant_output)
]


# Add EXACTLY the same map 100 times
for _ in range(100):

    constant_cmd.extend([
        "-cifti",
        str(constant_dlabel),
        "-index",
        "1"
    ])


print("\n========================================")
print("Creating constant 100-frame dlabel")
print("========================================")
print(f"Source: {constant_dlabel}")
print(f"Output: {constant_output}")

subprocess.run(constant_cmd, check=True)


print("\n========================================")
print("DONE")
print("========================================")

print("\nDynamic:")
print(dynamic_output)

print("\nConstant:")
print(constant_output)