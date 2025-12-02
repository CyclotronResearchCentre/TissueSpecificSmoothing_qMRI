#!/usr/bin/env python3
"""
qMRI Dataset: Aging Data from M.F. Callaghan
Equivalent Python script to MATLAB code for SUSAN bt list generation
"""

import os
import re
import glob
import numpy as np
import nibabel as nib
import pandas as pd

# --- Input dataset paths ---
# Adapt to your WSL path (no Windows drive letters!)
indir_bids_qMRI = "/home/antoi/smoothing/Data/BIDS_AgingData"
infol_qMRI = "AJ-TSPOON"

root_dir = os.path.join(indir_bids_qMRI, "derivatives", infol_qMRI)

# --- Read subject folders ---
sub_dirs = [d for d in os.listdir(root_dir) if d.startswith("sub-") and os.path.isdir(os.path.join(root_dir, d))]
nsub = len(sub_dirs)

# --- MPM maps to search ---
MPMs_listname = ["MTsat", "PDmap", "R1map", "R2starmap"]

wMPM_paths = [[] for _ in range(nsub)]
median_values = [[] for _ in range(nsub)]

for i, subj in enumerate(sub_dirs):
    subj_path = os.path.join(root_dir, subj, "anat")

    tmp_list = []
    tmp_values = []

    for modality in MPMs_listname:
        # Exclude gs_ files and match Mask_*MODALITY.nii
        regexp = re.compile(rf"^(?!.*gs_).*Mask_.*{modality}.*\.nii$")
        all_files = glob.glob(os.path.join(subj_path, "**", "*.nii"), recursive=True)

        file_list = [f for f in all_files if regexp.search(os.path.basename(f))]

        for nii_path in file_list:
            print(f"Sujet: {subj} | Modalité: {modality} | Fichier: {nii_path}")

            # --- Load NIfTI ---
            img = nib.load(nii_path).get_fdata()
            img_vec = img.flatten()

            # --- Count voxel categories ---
            n_pos = np.sum(img_vec > 0)
            n_zero = np.sum(img_vec == 0)
            n_neg = np.sum(img_vec < 0)
            n_nan = np.sum(np.isnan(img_vec))

            # --- Median ignoring zeros, negatives and NaNs (SUSAN convention) ---
            valid_vals = img_vec[(img_vec > 0) & ~np.isnan(img_vec)]
            med_val = np.median(valid_vals) * 0.75 if valid_vals.size > 0 else np.nan

            tmp_list.append(nii_path)
            tmp_values.append({
                "path": nii_path,
                "median075": med_val,
                "n_positive": int(n_pos),
                "n_zero": int(n_zero),
                "n_negative": int(n_neg),
                "n_nan": int(n_nan)
            })

    wMPM_paths[i] = tmp_list
    median_values[i] = tmp_values

print("Done. median_values now contains path → median×0.75 pairs.")

# --- Export to Excel/CSV ---
all_paths = []
all_median075 = []

for subj_vals in median_values:
    for entry in subj_vals:
        all_paths.append(entry["path"])
        all_median075.append(entry["median075"])

T = pd.DataFrame({"nii_path": all_paths, "median075": all_median075})

excel_path = os.path.join(root_dir, "SUSAN_bt_list.xlsx")
csv_path = os.path.join(root_dir, "SUSAN_bt_list.csv")

T.to_excel(excel_path, index=False)
T.to_csv(csv_path, index=False)

print(f"Export done: {excel_path} and {csv_path} created.")
