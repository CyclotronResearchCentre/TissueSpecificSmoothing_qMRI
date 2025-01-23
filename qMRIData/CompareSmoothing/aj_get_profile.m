function [profile, mni_coords] = aj_get_profile(matrix_paths, coord)
%--------------------------------------------------------------------------
% Function to extract a 1D signal at specific coordinates from a 3D matrix.
%
% INPUT
% matrix_paths:     List containing the paths to the 3D matrices where
%                   the profile extraction has to be done.
% coord:            Vector of 3 values including one NaN which specifies
%                   the coordinate along which the profile is taken.
%
% OUTPUT
% profile:          1D vector containing the values along the specified
%                   axis in the 3D matrix
% mni_coords:       Vector contaning the MNI coordinates of the extracted
%                   profile.
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
x = str2double(coord{1});
y = str2double(coord{2});
z = str2double(coord{3});

n = numel(matrix_paths);

profile = cell(1, n);
mni_coords = cell(1, n);  % Add associated MNI coordinates
for i = 1:n
    matrix = spm_read_vols(spm_vol(matrix_paths{i}));

    % Check which size is fixed
    if isnan(x) && ~isnan(y) && ~isnan(z)
        % x is variable: Extraction of the profile along x
        profile{i} = squeeze(matrix(:, y, z));
        mni_coords{i} = (1:size(matrix, 1)) - 90; % MNI coord for x
    elseif ~isnan(x) && isnan(y) && ~isnan(z)
        % y is variable: Extraction of the profile along y
        profile{i} = squeeze(matrix(x, :, z));
        mni_coords{i} = (1:size(matrix, 2)) - 126; % MNI coord for y
    elseif ~isnan(x) && ~isnan(y) && isnan(z)
        % z is variable: Extraction of the profile along z
        profile{i} = squeeze(matrix(x, y, :));
        mni_coords{i} = (1:size(matrix, 3)) - 72; % MNI coord for z
    else
        error('One and only one dimension should be left blank.');
    end
end

end
