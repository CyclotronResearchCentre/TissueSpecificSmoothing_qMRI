function final_img_info = aj_TS2one(x_paths, exMask_info, list_TC)
%--------------------------------------------------------------------------
% aj_TS2one - Combines tissue-specific images into one final image using 
% tissue-specific masks.
% This function takes tissue-weighted images (e.g., GM, WM, CSF) and 
% corresponding masks as inputs, processes them, and generates a final 
% combined image where each voxel value is selected based on the provided 
% masks. The function can be customized to handle any subset of tissues.
%
% INPUTS
% x_paths:      Cell array of file paths to the tissue-specific images 
%               (NIfTI files).
% exMask_info:  Cell array of NIfTI volume structures for tissue masks.
% list_TC:      (Optional) Number or list of tissue classes of interest. 
%               If not provided, all available tissue classes will be used.
%
% OUTPUT
% final_img_info: Structure of the final combined image, saved as a NIfTI
%               file.
%
% PROCESS
% 1. Validates input paths and loads image data into `x_info`.
% 2. Checks that all images have matching dimensions.
% 3. Validates and resizes tissue masks to match the image dimensions if 
%    needed.
% 4. Reads the voxel data from tissue images and masks.
% 5. Applies mask conditions to build a final image with appropriate voxel 
%    values.
% 6. Writes the final combined image to a NIfTI file.
%
% LIMITATIONS
% - The function requires that all input images and masks have the same
%   spatial dimensions and resolution. Resizing is performed with nearest- 
%   neighbor interpolation, which may not be ideal for all applications.
% - The function assumes binary masks (values 0 or 1) and applies strict 
%   conditions for combining images.
% - The process may be limited by available memory for large images.
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------

% Dealing with inputs
if nargin < 2
    error('aj_TS2one ISSUE: not enough inputs.\n');
end

% Load x_paths into x_info and validate
if iscell(x_paths)
    nb_img = length(x_paths); % number of images of different tissue class
    x_info = cell(1,nb_img);
    for i = 1:nb_img
        % Check if x_paths is a valid file path
        if exist(x_paths{i}, 'file') == 2  % '2' indicates a file
            x_info{i} = spm_vol(x_paths{i});
        else
            error('aj_TS2one ISSUE: input is not an existing nifti file.');
        end

        if ~isstruct(x_info{i}) || ~all(isfield(x_info{i}, {'fname', 'dim', 'mat'}))
            fprintf('aj_TS2one ISSUE: invalid structure %d in cell array x_info.\n', i);
            return;
        end
    end
else
    warning('aj_TS2one ISSUE: x_info is not a cell array.\n');
    return;
end

% Check if all images have the same dimensions
for i = 1:nb_img
    for j = 1:nb_img
        if i ~= j && ~isequal(x_info{i}.dim, x_info{j}.dim)
            fprintf('aj_TS2one ISSUE: not same matrix dimensions for files %s and %s.\n', x_info{i}.fname, x_info{j}.fname);
            return;
        end
    end
end

% Load and validate exMask_info
if iscell(exMask_info)
    nb_mask = length(exMask_info); % number of masks of different tissue class
    for i = 1:nb_mask
        if ~isstruct(exMask_info{i}) || ~all(isfield(exMask_info{i}, {'fname', 'dim', 'mat'}))
            fprintf('aj_TS2one ISSUE: invalid structure %d in cell array exMask_info.\n', i);
            return;
        end
        
        % Check if mask dimensions match
        for j = 1:nb_img
            if i ~= j && ~isequal(exMask_info{i}.dim, exMask_info{j}.dim) % Compare if the current file dimensions match with others
                fprintf('aj_TS2one ISSUE: not same matrix dimensions for files %s and %s.\n', exMask_info{i}.fname, exMask_info{j}.fname);
                return;
            end
        end
    end
else
    warning('aj_TS2one ISSUE: exMask_info is not a cell array.\n');
    return;
end

% Use all tissue classes if no specific list is provided
if nargin < 3
    fprintf('aj_TS2one: no provided list of tissue classes of interest (list_TC). All tissue classes will be considered (default option).\n');
    list_TC = nb_img;
end

if nb_mask < list_TC
    error('aj_TS2one ISSUE: nb_mask < list_TC.');
elseif nb_img < list_TC
    error('aj_TS2one ISSUE: nb_img < list_TC.');
else
    nb_TC = list_TC;
end

%--------------------------------------------------------------------------
% Check if masks and x_img have the same dimensions and build a mask 
% having the same dimensions than x_img if needed
dim_exMask_info = cell(1, nb_TC);

if ~isequal(x_info{1}.dim, exMask_info{1}.dim)
    fprintf('Resizing explicit masks to good dimensions...\n');
    
    % Create a unity matrix of the same dimensions as x_info
    unity_matrix = ones(x_info{1}.dim);
    dim_exMask_info{i} = spm_file(exMask_info{i}.fname, 'prefix', 'unity_');
    temp_vol = exMask_info{i}; % Use an existing volume for header template
    temp_vol.fname = dim_exMask_info{i}; % Set the output file name
    temp_vol.dim = x_info{1}.dim; % Adjust dimensions to match x_info
    temp_vol.mat = x_info{1}.mat; % Ensure alignment of the coordinate system
    unity_vol = spm_write_vol(temp_vol, unity_matrix); % % Write unity matrix to a NIfTI file
    
    % Apply spm_imcalc to resize masks
    ic_flag = struct('interp', 0); % Nearest-neighbor interpolation for binary masks
    for i = 1:nb_TC
        Vi = [unity_vol, exMask_info{i}];
        Vo_in = spm_file(exMask_info{i}.fname, 'prefix', 'resized_');
        f = '(i1./i2)';  
        dim_exMask_info{i} = spm_imcalc(Vi, Vo_in, f, ic_flag);
        
        if isequal(dim_exMask_info{i}.dim, x_info{1}.dim)
            fprintf('aj_TS2one: well resized mask in %s.\n', dim_exMask_info{i}.fname);
        else
            fprintf('aj_TS2one ISSUE: still not the same dimension for %s.\n', dim_exMask_info{i}.fname);
            return;
        end
    end
    
    delete(unity_vol.fname);
    
else
    dim_exMask_info = exMask_info;
end

% Load the explicit masks
dim_exMask = cell(1,nb_TC);
for i = 1:nb_TC
    dim_exMask{i} = spm_read_vols(dim_exMask_info{i});
end

final_image = zeros(size(dim_exMask{1}));  % Initialize the final image

% Loop through each specified tissue class (GM, WM, etc.)
for i = 1:nb_TC  % Loop over the selected tissue-weighted images
    % Read the NIfTI file of the current tissue image
    x = spm_read_vols(spm_vol(x_info{i}));  
    
    % Create a mask condition for the current tissue class
    mask_condition = dim_exMask{i} == 1;  % The current mask must be 1 (active)

    % Ensure that all other masks are not active (set to 0)
    for j = 1:nb_TC
        if j ~= i
            mask_condition = mask_condition & (dim_exMask{j} == 0);
        end
    end
    
    % Apply the mask condition to set the voxel values in the final image
    final_image(mask_condition) = x(mask_condition);
end

% Save the final image to a NIfTI file
final_img_info = x_info{1};  % Use the header from the first image (x1) as reference
final_img_info.fname = fullfile(fileparts(x_paths{1}), 'final_result.nii');
final_img_info.descrip = 'Final result based on explicit masks (GM, WM, CSF)';
spm_write_vol(final_img_info, final_image);

fprintf('Final image has been created and saved as %s\n', final_img_info.fname);

end
