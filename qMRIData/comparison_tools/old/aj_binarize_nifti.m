%% Binarize NIfTI Image
%
% -------------------------------------------------------------------------
% DESCRIPTION
% -------------------------------------------------------------------------
% This function converts a NIfTI image into a binary mask.
%
% All voxels satisfying:
%       voxel ~= 0 AND voxel is not NaN
% are assigned a value of 1.
%
% All remaining voxels are assigned a value of 0.
%
% The resulting binary image is saved as a new NIfTI file using SPM.
%
% This helper function is particularly useful for:
%   - Thresholded statistical maps
%   - Significant voxel masks
%   - SnPM/SPM FWE-corrected images
%   - Cluster-based analyses
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
% infile : char
%     Full path to the input NIfTI image.
%
% outfile : char
%     Full path to the output binary NIfTI image.
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
% A binary NIfTI image is written to disk:
%
%     1 = significant/non-zero valid voxel
%     0 = background, zero or NaN voxel
%
% The output datatype is uint8-compatible integer format (SPM datatype 2).
%
% -------------------------------------------------------------------------
% REQUIREMENTS
% -------------------------------------------------------------------------
% - SPM12 must be added to the MATLAB path.
%
% -------------------------------------------------------------------------
% AUTHOR
% -------------------------------------------------------------------------
% Antoine Jacquemin
%
% -------------------------------------------------------------------------
% EXAMPLE
% -------------------------------------------------------------------------
% infile  = 'Tmap_thr.nii';
% outfile = 'Tmap_thr_bin.nii';
%
% aj_binarize_nifti(infile, outfile);
%
% -------------------------------------------------------------------------

function aj_binarize_nifti(infile, outfile)

    %----------------------------------------------------------------------
    % Load input NIfTI image
    %----------------------------------------------------------------------
    V = spm_vol(infile);
    Y = spm_read_vols(V);

    %----------------------------------------------------------------------
    % Binarization
    %----------------------------------------------------------------------
    % Keep only non-zero and non-NaN voxels
    Ybin = double((Y ~= 0) & ~isnan(Y));

    %----------------------------------------------------------------------
    % Define output image properties
    %----------------------------------------------------------------------
    V.fname = outfile;

    % SPM datatype = 2 → unsigned 8-bit integer
    V.dt = [2 0];

    %----------------------------------------------------------------------
    % Save binary image
    %----------------------------------------------------------------------
    spm_write_vol(V, Ybin);

    fprintf('Binary image saved: %s\n', outfile);

end