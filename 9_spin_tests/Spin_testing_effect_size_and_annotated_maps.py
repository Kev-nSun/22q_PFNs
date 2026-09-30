from pathlib import Path

import nibabel as nib
import pandas as pd
import numpy as np

from neuromaps import datasets, nulls, transforms
from neuromaps.stats import compare_images

# Spin tests require vertex counts to be in 32492 per hemisphere

# ----------------------------
# Settings
# ----------------------------
N_PERM = 5000
SEED = 1234
METRIC = "pearsonr"

results_dir = Path(
    "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/"
    "Analyses/Results/Spin_Tests/"
)
results_dir.mkdir(parents=True, exist_ok=True)

nulls_dir = results_dir / "null_distributions"
nulls_dir.mkdir(parents=True, exist_ok=True)


# ----------------------------
# Read effect-size and variability map
# ----------------------------
TCA_area_indir = Path(
    "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/"
    "Analyses/Results/TCA_vert_SA_gam_results/Maps/Unthresholded/GIFTI/"
)
TCA_area_map = nib.load( # type: ignore
    str(TCA_area_indir / "Rel_1_logTCA_vert_SMOOTH_logSA_gam_weights_unthresholded.gii")
)

TCA_PFNs_indir = Path(
    "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/"
    "Analyses/Results/TCA_PFN_loadings_gam_results/NORMED/Maps/Unthresholded/to-be-GIFTIed/GIFTI/"
)
TCA_PFN_map = nib.load( # type: ignore
    str(TCA_PFNs_indir / "ABS_MAX_logTC_PFN_loading_gam_weights_LOG1.5_unthres.gii")
)

Del_area_indir = Path(
    "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/"
    "Analyses/Results/22q_vert_SA_gam_results/Maps/Unthresholded/Final_map/GIFTI/"
)
Del_area_map = nib.load( # type: ignore
    str(Del_area_indir / "22q_vert_SA_SMOOTH_subclass_gam_weights_EXP_REL_1_unthresholded.gii")
)

Del_PFNs_indir = Path(
    "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/"
    "Analyses/Results/22q_PFN_loadings_gam_results_subclass/Normed_TCA/Maps/Unthresholded/to-be-GIFTIed/GIFTI/"
)
Del_PFN_map = nib.load( # type: ignore
    str(Del_PFNs_indir / "ABS_MAX_22q_PFN_loading_gam_weights_unthres.gii")
)

MAD_indir = Path(
    "C:/Users/kevin/OneDrive/Documents/NGG_PhD/Alexander-Bloch/22q_Project/"
    "Analyses/Results/Spin_Tests/MAD_&_resid_maps/GIFTI"
)
TCA_PFN_resid_map = nib.load( # type: ignore
    str(MAD_indir / "Scaling_MAD_resid_normed.gii")
)
Del_PFN_resid_map = nib.load( # type: ignore
    str(MAD_indir / "del_22q_MAD_resid_normed.gii")
)
MAD_map = nib.load( # type: ignore
    str(MAD_indir / "MAD_Avg.gii")
)

effect_maps = {
    "TCA_area": TCA_area_map,
    "Del_area": Del_area_map,
    "TCA_PFN": TCA_PFN_map,
    "Del_PFN": Del_PFN_map,
    "TCA_PFN_resid": TCA_PFN_resid_map,
    "Del_PFN_resid": Del_PFN_resid_map,
    "MAD_var": MAD_map,
}


# ----------------------------
# Fetch / transform annotation maps to fsLR 32k
# ----------------------------
princ_grad = datasets.fetch_annotation(
    source="margulies2016",
    desc="fcgradient01",
    space="fsLR",
    den="32k",
)

# FC_grad2 = datasets.fetch_annotation(
#     source="margulies2016",
#     desc="fcgradient02",
#     space="fsLR",
#     den="32k",
# )

abagen_PC1 = datasets.fetch_annotation(
    source="abagen",
    desc="genepc1",
    space="fsaverage",
    den="10k",
)
abagen_PC1_fslr = transforms.fsaverage_to_fslr(abagen_PC1, "32k")

# dev_exp = datasets.fetch_annotation(
#     source="hill2010",
#     desc="devexp",
#     space="fsLR",
#     den="164k",
# )
# dev_exp_32k = transforms.fslr_to_fslr(dev_exp, "32k")

# evo_exp = datasets.fetch_annotation(
#     source="hill2010",
#     desc="evoexp",
#     space="fsLR",
#     den="164k",
# )
# evo_exp_32k = transforms.fslr_to_fslr(evo_exp, "32k")

myelin = datasets.fetch_annotation(
    source="hcps1200",
    desc="myelinmap",
    space="fsLR",
    den="32k",
)

cort_thickness = datasets.fetch_annotation(
    source="hcps1200",
    desc="thickness",
    space="fsLR",
    den="32k",
)

neurosynth_PC1 = datasets.fetch_annotation(
    source="neurosynth",
    desc="cogpc1",
    space="MNI152",
    res="2mm",
)
neurosynth_PC1_fslr = transforms.mni152_to_fslr(neurosynth_PC1, "32k")

gluc_met = datasets.fetch_annotation(
    source="raichle",
    desc="cmrglc",
    space="fsLR",
    den="164k",
)
gluc_met_32k = transforms.fslr_to_fslr(gluc_met, "32k")

CBF = datasets.fetch_annotation(
    source="satterthwaite2014",
    desc="meancbf",
    space="MNI152",
    res="1mm",
)
CBF_fslr = transforms.mni152_to_fslr(CBF, "32k")

# fc_homo = datasets.fetch_annotation(
#     source="xu2020",
#     desc="FChomology",
#     space="fsLR",
#     den="32k",
# )

evo_exp_xu = datasets.fetch_annotation(
    source="xu2020",
    desc="evoexp",
    space="fsLR",
    den="32k",
)

annot_maps = {
    "princ_grad": princ_grad,
    # "FC_grad2": FC_grad2,
    "abagen_PC1": abagen_PC1_fslr,
    # "dev_exp": dev_exp_32k,
    # "evo_exp": evo_exp_32k,
    "myelin": myelin,
    "cort_thickness": cort_thickness,
    "neurosynth_PC1": neurosynth_PC1_fslr,
    "gluc_met": gluc_met_32k,
    "CBF": CBF_fslr,
    # "fc_homo": fc_homo,
    "evo_exp_xu": evo_exp_xu,
}


# ----------------------------
# Select maps that will be spun
# ----------------------------
spin_source_names = {
    "TCA_area",
    "Del_area",
    "TCA_PFN",
    "Del_PFN",
    "TCA_PFN_resid",
    "Del_PFN_resid",
    "MAD_var",
}

spin_maps = {
    name: effect_maps[name]
    for name in spin_source_names
}

# Optional validation
missing_sources = spin_source_names.difference(effect_maps)

if missing_sources:
    raise KeyError(
        f"These spin source maps are missing from effect_maps: "
        f"{sorted(missing_sources)}"
    )


# ----------------------------
# Precompute nulls only for selected source maps
# ----------------------------
spins = nulls.alexander_bloch(
    None,
    atlas="fsLR",
    density="32k",
    n_perm=N_PERM,
    seed=SEED,
)

effect_map_nulls = {}

for src_name, src_map in spin_maps.items():

    effect_map_nulls[src_name] = nulls.alexander_bloch(
        src_map,
        atlas="fsLR",
        density="32k",
        spins=spins,
    )


# ----------------------------
# Run comparisons
# ----------------------------
for src_name, src_map in spin_maps.items():
    rows = []
    src_nulls = effect_map_nulls[src_name]

    # Compare spun source to every other effect-size map
    for trg_name, trg_map in effect_maps.items():
        if trg_name == src_name:
            continue

        similarity, pvalue, null_dist = compare_images( # type: ignore
            src_map,
            trg_map,
            metric=METRIC,
            nulls=src_nulls,
            return_nulls=True,
        )

        null_out = np.empty((len(null_dist) + 1, 2), dtype=object)
        null_out[0] = ["observed", similarity]
        null_out[1:, 0] = "null"
        null_out[1:, 1] = null_dist

        np.savetxt(
            nulls_dir / f"{src_name}_vs_{trg_name}_nulls.csv",
            null_out,
            delimiter=",",
            fmt="%s",
        )

        rows.append({
            "source_map": src_name,
            "target_map": trg_name,
            "target_type": "effect_map",
            "metric": METRIC,
            "n_perm": N_PERM,
            "r": float(similarity),
            "p_spin": float(pvalue),
        })

    # Compare spun source to every annotation map
    for trg_name, trg_map in annot_maps.items():
        similarity, pvalue, null_dist = compare_images( # type: ignore
            src_map,
            trg_map,
            metric=METRIC,
            nulls=src_nulls,
            return_nulls=True,
        )

        null_out = np.empty((len(null_dist) + 1, 2), dtype=object)
        null_out[0] = ["observed", similarity]
        null_out[1:, 0] = "null"
        null_out[1:, 1] = null_dist

        np.savetxt(
            nulls_dir / f"{src_name}_vs_{trg_name}_nulls.csv",
            null_out,
            delimiter=",",
            fmt="%s",
        )

        rows.append({
            "source_map": src_name,
            "target_map": trg_name,
            "target_type": "annotation",
            "metric": METRIC,
            "n_perm": N_PERM,
            "r": float(similarity),
            "p_spin": float(pvalue),
        })

    out_df = pd.DataFrame(rows)

    out_df.to_csv(
        results_dir / f"{src_name}_spin_results_w_MAD_5000.csv",
        index=False,
    )

print("Done.")