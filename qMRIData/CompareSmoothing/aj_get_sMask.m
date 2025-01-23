function sMasks_paths = aj_get_sMask(gs_TCmap_paths, pth_out)
%--------------------------------------------------------------------------
% Function to rebuild the smoothed masks which are used as TSPOON
% denominator from the TWS denominator.
%
% INPUT
% gs_TCmap_paths:   List containing the paths to the smoothed
%                   tissue-specific maps (i.e. TWS denominator)
% pth_out:          Output path to save the results.
%
% OUTPUT
% sMasks_paths:     List of paths to the smoothed masks
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
fwhm = 3; % ref to Callaghan article
nTC = numel(gs_TCmap_paths);

% Get the explicit tissue-specific binary masks from gs_mwc-images
% -> Mask_gs_mwc-images
addpath('C:\Users\antoi\Documents\master_thesis\MATLAB\smoothing');
Mask_gs_mwc_paths = aj_exMask(gs_TCmap_paths, pth_out);

% Gaussian isotropic smoothing of Mask_gs_mwc-images
% -> gs_Mask_gs_mwc-images
gs_Mask_gs_mwc_paths = cell(1,nTC);
for jj = 1:nTC
    gs_Mask_gs_mwc_paths{jj} = spm_file(Mask_gs_mwc_paths{jj},'prefix','gs_');
    spm_smooth(Mask_gs_mwc_paths{jj}, gs_Mask_gs_mwc_paths{jj}, fwhm);
end

% Managing small values of gs_Mask_gs_mwc-images to avoid dividing by 0 or
% too small values by both creating a binary mask at threshold of 5% (non 
% linear operation) and converting 0 to NaN.
% -> thr_gs_Mask_gs_mwc-images = d-images
thr_gs_Mask_gs_mwc_paths = cell(1,nTC);
for jj = 1:nTC
    Vi = char(gs_Mask_gs_mwc_paths{jj});
    thr_gs_Mask_gs_mwc_paths{jj} = spm_file(gs_Mask_gs_mwc_paths{jj},'prefix','thr_','path',pth_out);
    f = 'i1 .* (i1>.05)';
    % Flags for image calculation
    ic_flag = struct(... % type is set below based on that of input image
        'interp', -4);   % 4th order spline interpolation
    spm_imcalc(Vi, thr_gs_Mask_gs_mwc_paths{jj}, f, ic_flag);

    % Set 0's to NaN for the types that have NaN representation
    % -> proper implicit masking during SPM analysis
    if spm_type(spm_vol(gs_TCmap_paths{1}).dt(1),'nanrep') % check if mwc* has NaN representation
        hmri_proc_zero2nan(thr_gs_Mask_gs_mwc_paths{jj});
    end
end

sMasks_paths = cell(1,3);
sMasks_paths{1} = Mask_gs_mwc_paths;
sMasks_paths{2} = gs_Mask_gs_mwc_paths;
sMasks_paths{3} = thr_gs_Mask_gs_mwc_paths;

fprintf('All TSPOON masks at each step are successfully saved in %s.\n', pth_out);

end