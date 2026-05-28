function mergedImg_path = aj_TS2One(x_paths, exMask_paths, keeping_names)
%--------------------------------------------------------------------------
% aj_TS2One - Combines tissue-specific images into one final image using 
% tissue-specific masks.
% This function takes tissue-weighted images (e.g., GM, WM, CSF) and 
% corresponding masks as inputs, processes them, and generates a final 
% combined image where each voxel value is selected based on the provided 
% masks. The function can be customized to handle any subset of tissues.
%
% INPUT
% x_paths:      Cell array of file paths to the tissue-specific images 
%               (NIfTI files).
% exMask_info:  Cell array of NIfTI file paths to the tissue masks.
%
% OUTPUT
% mergedImg_path: Path to the merged tissue-classes image, saved as a NIfTI
%               file.
%
% LIMITATIONS
% - The function requires that all input images and masks have the same
%   spatial dimensions and resolution. Resizing is performed with nearest- 
%   neighbor interpolation, which may not be ideal for all applications.
% - The function assumes binary masks (values 0 or 1) and applies strict 
%   conditions for combining images.
% - Input files have to be sorted (GM, WM etc)
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
%% Dealing with inputs
if nargin < 2
    error('aj_TS2One ISSUE: not enough inputs.\n');
end

if iscell(x_paths) || iscell(exMask_paths)
    nx = numel(x_paths);
    nMask = numel(exMask_paths);
    if nx<2 || nMask<2
        error('At least one input contains only one cell for one tissue-specific maps.');
    elseif nx ~= nMask
        error('Not the same number of masks and tissue-specific masks.');
    end
else
    error('At least one input is not a cell or cell array.');
end

%% Do the job
% Flags for image calculation
ic_flag = struct('interp', -4);

Vi = char(x_paths{1},exMask_paths{1});
maskGM_path = spm_file(x_paths{1},'prefix','GM_');
f = '(i1.*i2)';
spm_imcalc(Vi, maskGM_path, f, ic_flag);

Vi = char(x_paths{2},exMask_paths{2});
maskWM_path = spm_file(x_paths{2},'prefix','WM_');
f = '(i1.*i2)';
spm_imcalc(Vi, maskWM_path, f, ic_flag);

Vi = char(maskGM_path, maskWM_path);

current_name = spm_file(x_paths{1}, 'basename');
nNames = numel(keeping_names);
for i = 1:nNames
    result = contains(current_name, keeping_names{i});
    if result
        new_name = sprintf('merged_%s.nii',keeping_names{i});
    end
end

current_fullPath = spm_file(x_paths{1},'fpath');
mergedImg_path = fullfile(current_fullPath, new_name);
f = '(i1 + i2)';
spm_imcalc(Vi, mergedImg_path, f, ic_flag);

[pth,nam,ext,~] = spm_fileparts(maskGM_path);
delete(fullfile(pth,[nam,ext]));
[pth,nam,ext,~] = spm_fileparts(maskWM_path);
delete(fullfile(pth,[nam,ext]));

fprintf('Final image has been created and saved as %s\n', mergedImg_path);

end