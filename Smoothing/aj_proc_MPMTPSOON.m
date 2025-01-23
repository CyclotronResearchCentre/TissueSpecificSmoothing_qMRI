function [fn_out,fn_sMask] = aj_proc_MPMTPSOON(fn_wMPM, fn_mwTC, fwhm, l_TC, pth_out)
%--------------------------------------------------------------------------
% This function applies T-SPOON (Tissue-SPecific smOOthing cOmpensated) to 
% multiple parametric and tissue class maps from a single subject. It 
% smooths modulated, warped tissue class maps, generates tissue-specific 
% binary masks, applies a threshold to avoid divisions by zero, and finally 
% outputs the tissue-specific smoothed parametric maps (T-SPOON images).
%
% FORMAT
% [fn_out, fn_smwc] = aj_proc_MPMTPSOON(fn_wMPM, fn_mwTC, fwhm, l_TC, pth_out)
% 
% INPUTS
% - fn_wMPM : filenames (char array) of the warped MPM, i.e. the
%             w*MT/R1/R2s/A.nii files
% - fn_mwTC : filenames (char array) of the modulated warped tissue
%             classes, i.e. the mwc1/2*.nii files
% - fwhm    : width of smoothing kernel in mm [6 by def.]
% - l_TC    : explicit list of tissue classes used [1:nTC by def.]
% - pth_out : output path [empty by default -> same as input maps]
% 
% OUTPUTS
% - fn_out  :   Cell array (one cell per MPM) of filenames (char array) 
%               containing the "TSPOON smoothed tissue specific MPMs".
% - fn_sMask:   Char array of paths of (non thresholded) smoothed masks of   
%               modulated warped tissue classes.
%
% PROCESS:
% 1. Checks and sets input parameters and paths.
% 2. Smooths the input tissue class maps using Gaussian isotropic smoothing.
% 3. Generates tissue-specific binary masks and applies a threshold.
% 4. Computes tissue-specific smoothed parametric maps by dividing smoothed
%    parametric images by corresponding masks.
% 5. Ensures that zero values in the masks are converted to NaN to prevent 
%    division errors and maintain valid data.
% 6. Outputs the T-SPOON images and deletes temporary intermediate files.
%
% REFERENCES
% J. E. Lee & al. (2008), https://doi.org/10.1016/j.neuroimage.2008.09.041
%--------------------------------------------------------------------------
% NOTES ABOUT DENOMINATOR (d-images)
% d-images = the thresholded, Gaussian isotropic smoothed, tissue-specific
% binary masks of Gaussian isotropic smoothed, subject-specific modulated 
% warped tissue class maps. "Subject-specific" so it can be computed once
% and used for all MPMs of a same subject. d-images are computed for all
% available tissue classes in fn_mwTC. The tissue-specific binary masks are
% computed by follwing the "majority and upper than 20% threshold".
% Majority implies taking into account all the available tissue classes in
% fn_mwTC.
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Deal with inputs
% Check input
if nargin<5, pth_out = []; end
if nargin<3, fwhm = 6; end
if nargin<2
    error('hMRI: TSPOON smoothing: Provide 2 inputs, see help.');
end

% Count images and check
nMPM = size(fn_wMPM,1);
nTC  = size(fn_mwTC,1);

TC_names = {'GM', 'WM', 'CSF'};
nTCNames = length(TC_names);

if nTCNames<nTC
    error('hMRI: TSPOON smoothing: List of TC names not long enough.');
end

if nargin<4 % changed from TWS
    l_TC = 1:nTC;
end

% change from TWS
nl_TC = max(l_TC);

% Flags for image calculation
ic_flag = struct(... % type is set below based on that of input image
    'interp', -4);   % 4th order spline interpolation

%% Computing d-images
% ic_flag for d-images computes for the first fn_wMPM. However d-images are
% used for all fn_wMPM (doen't matter in TWS because no use of spm_imcalc).
% Hypothesis:
% 1. not important for d-images because correct if needed during the
% division n/d : FALSE (wrong results)
% 2. test for to be fn_mwTC-specific
V_ii = spm_vol(fn_wMPM(1,:));
ic_flag.dtype = V_ii.dt(1);

% Gaussian isotropic smoothing of modulated warped tissue class maps 
% -> gs_mwc-images
gs_mwc_paths = cell(1,nTC);
for jj=1:nTC
    gs_mwc_paths{jj} = spm_file(fn_mwTC(jj,:),'prefix','gs_','number','','path',pth_out); %even if func and not anat
    spm_smooth(fn_mwTC(jj,:),gs_mwc_paths{jj},fwhm); % smooth mwc(jj)
end

% Get the explicit tissue-specific binary masks from gs_mwc-images
% -> Mask_gs_mwc-images
Mask_gs_mwc_paths = aj_exMask(gs_mwc_paths, pth_out);

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
    spm_imcalc(Vi, thr_gs_Mask_gs_mwc_paths{jj}, f, ic_flag);

    % Set 0's to NaN for the types that have NaN representation
    % -> proper implicit masking during SPM analysis
    if spm_type(spm_vol(fn_mwTC(jj,:)).dt(1),'nanrep') % check if mwc* has NaN representation
        hmri_proc_zero2nan(thr_gs_Mask_gs_mwc_paths{jj});
    end
end

%% Do the job
% Initialize output and loop over MPMs
fn_out = cell(nMPM,1);
for ii=1:nMPM
    % ii^th MPM to be treated
    fn_wMPM_ii = fn_wMPM(ii,:);
    
    % Deal with image type -> keep original
    % + find out if it has NaN representation
    V_ii = spm_vol(fn_wMPM_ii);
    ic_flag.dtype = V_ii.dt(1);
    NaNrep = spm_type(V_ii.dt(1),'nanrep');
        
    % Get TC-segmented/masked MPM -> MaskMPM-images
    MaskMPM_paths = cell(1,nl_TC);
    for jj = 1:nl_TC
        Vi = char(fn_wMPM_ii, Mask_gs_mwc_paths{jj}); % spm_imcalc will keep fn_wMPM dimensions (first input)
        MaskMPM_paths{jj} = spm_file(fn_wMPM_ii,'prefix',['Mask',num2str(jj),'_'],'path',pth_out);
        f = '(i1.*i2)';
        spm_imcalc(Vi, MaskMPM_paths{jj}, f, ic_flag);
    end
    
    % Gaussian isotropic smoothing of the MaskMPMs -> gs_MaskMPM-images
    gs_MaskMPM_paths = cell(1,nl_TC);
    for jj = 1:nl_TC
        gs_MaskMPM_paths{jj} = spm_file(MaskMPM_paths{jj},'prefix','gs_');
        spm_smooth(MaskMPM_paths{jj}, gs_MaskMPM_paths{jj}, fwhm);
    end
    
    % Dividing gs_MaskMPMs by xxx to get T-SPOON signals -> tspoon-images
    tspoon_paths = cell(1,nl_TC);
    for jj = 1:nl_TC
        Vi = char(gs_MaskMPM_paths{jj}, thr_gs_Mask_gs_mwc_paths{jj}); % spm_imcalc will keep fn_wMPM dimensions
        tspoon_paths{jj} = spm_file(fn_wMPM_ii,'prefix',['TSPOON_',TC_names{jj},'_'],'path',pth_out);
        f = '(i1./i2)';
        spm_imcalc(Vi, tspoon_paths{jj}, f, ic_flag);
        
        % Set 0's to NaN for the types that have NaN representation
        % -> proper implicit masking during SPM analysis
        if NaNrep
            hmri_proc_zero2nan(tspoon_paths{jj});
        end
    end
    
    fn_out{ii} = char(tspoon_paths); % saved as char array
  
    % Delete unuseful files
    fn_2delete = char(char(MaskMPM_paths),char(gs_MaskMPM_paths));
    for jj=1:size(fn_2delete,1)
        delete(deblank(fn_2delete(jj,:)));
    end
end

fn_sMask = char(gs_Mask_gs_mwc_paths); % catch the non thresholded smoothed masks of smwTC images

% Delete unuseful files
fn_2delete = char(char(Mask_gs_mwc_paths), char(gs_mwc_paths), char(thr_gs_Mask_gs_mwc_paths));
for jj=1:size(fn_2delete,1)
    delete(deblank(fn_2delete(jj,:)));
end
end