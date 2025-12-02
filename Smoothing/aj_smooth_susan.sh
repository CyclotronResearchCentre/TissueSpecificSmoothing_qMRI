/////////////////////////////////////////////////////////////////////////

nano ~/.bashrc
export PATH=$PATH:/home/antoi/smoothing/FSL-susan
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/antoi/smoothing/FSL-susan/znzlib
source ~/.bashrc
which susan
conda activate smoothing

export FSLOUTPUTTYPE=NIFTI_GZ

/////////////////////////////////////////////////////////////////////////

#!/bin/bash
# Bash script to run SUSAN for each NIfTI listed in CSV

CSV="/home/antoi/smoothing/Data/BIDS_AgingData/derivatives/AJ-TSPOON/SUSAN_bt_R1.csv"
FWHM=3.0

# Loop through CSV, skipping header
tail -n +2 "$CSV" | while IFS=, read -r nii_path median075; do
    # Check file exists
    if [[ ! -f "$nii_path" ]]; then
        echo "❌ File not found: $nii_path"
        continue
    fi

    # Compute dt from FWHM (sigma = gaussian standard deviation)
    dt=$(python3 -c "import numpy as np; print($FWHM/np.sqrt(8*np.log(2)))")

    out_nii="$(dirname "$nii_path")/susan_$(basename "${nii_path%.nii.gz}" .nii).nii.gz"

    echo "--------------------------------------------------"
    echo "Input:  $nii_path"
    echo "BT:     $median075"
    echo "DT:     $dt"
    echo "Output: $out_nii"
    echo "--------------------------------------------------"

    # Run SUSAN
    susan "$nii_path" "$median075" "$dt" 3 0 0 "$out_nii"
done

/////////////////////////////////////////////////////////////////////////

susan --help

gzip -c /home/antoi/smoothing/Data/BIDS_AgingData/derivatives/AJ-TSPOON/sub-001/anat/Mask_GM_sub-S001_space-MNI_MTsat.nii \
  > /home/antoi/smoothing/Data/BIDS_AgingData/derivatives/AJ-TSPOON/sub-001/anat/Mask_GM_sub-S001_space-MNI_MTsat.nii.gz
