function [gs_imgaussfilt3_paths, gs_spm_paths,...
    tws_paths, smwTC_paths,...
    tspoon_paths, sMask_paths] = ...
    aj_smoothing(MPM_paths, TCseg_paths, param, flag, dim)
%--------------------------------------------------------------------------
% Function to apply different smoothing techniques (Gaussian, TWS, TSPOON)
% to the given data, based on the flags provided.
%
% INPUTS
% MPM_paths:    A list containing paths to images to smooth
% TCseg_paths:  A list containing paths to tissue-specific segmented images
% param:        Smoothing parameters (from default)
% flag:         Flags to determine which smoothing techniques to apply
%               (from default)
% dim:          Number for dimension of the data (1 for 1D, 2 for 2D or 3
%               for 3D) 
%
% OUTPUTS
% Default output are 0 (depending on flags)
% gsP_signal:       Gaussian smoothed signal from imgaussfilt3
% gs_path:          Gaussian smoothed file path from spm_smooth
% twsP_signal:      Tissue-weighted smoothed signal
% smwTC_paths:      Paths to the smoothed modulated warped tissue class
%                   maps
% tspoon_paths:     Cell array of filenames (nifti files) of the Tissue-
%                   SPecific smOOthing compeNsated (TSPOON) for GM & WM
% sMask_paths:      Cell array of filenames (nifti files) of the individual
%                   explicit mask
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
gs_imgaussfilt3_paths = 0;
gs_spm_paths = 0;
tws_paths = 0;
smwTC_paths = 0;
tspoon_paths = 0;
sMask_paths = 0;

%% Apply standard Gaussian smoothing
if flag.gaussian
    disp('Executing standard Gaussian smoothing...');
    
    % Automatically create a new smoothing-specific derivatives directory for the subject
    pth_out = regexprep(spm_file(MPM_paths(1,:), 'fpath'), param.outDerivName, 'AJ-GS');
    if ~exist(pth_out, 'dir')
        mkdir(pth_out);
    else
        % If the directory already exists, delete all the files inside it
        file_list = dir(fullfile(pth_out, '*'));
        for k = 1:length(file_list)
            % Ignore the special directories '.' and '..'
            if ~file_list(k).isdir
                delete(fullfile(pth_out, file_list(k).name));
            end
        end
    end
    
    [gs_imgaussfilt3_paths, gs_spm_paths] = aj_smooth_gaussian(MPM_paths, pth_out, param.fwhm_gs, dim);
end

%% Apply tissue-weighted smoothing (TWS)
if flag.tws
    disp('Executing tissue-weighted smoothing (TWS)...');
    
    % Taking the tissue probability map paths from TPM.nii
    nTC = size(TCseg_paths,1);
    TPM_paths = cell(1,nTC);
    for i = 1:nTC
        TPM_paths{i} = fullfile(spm('Dir'), 'tpm', ['TPM.nii,' num2str(i)]);
    end
    
    % Automatically create a new smoothing-specific derivatives directory for the subject
    pth_out = regexprep(spm_file(MPM_paths(1,:), 'fpath'), param.outDerivName, 'AJ-TWS');
    if ~exist(pth_out, 'dir')
        mkdir(pth_out);
    else
        % If the directory already exists, delete all the files inside it
        file_list = dir(fullfile(pth_out, '*'));
        for k = 1:length(file_list)
            % Ignore the special directories '.' and '..'
            if ~file_list(k).isdir
                delete(fullfile(pth_out, file_list(k).name));
            end
        end
    end
    
    [tws_paths, smwTC_paths] = hmri_proc_MPMsmooth(char(MPM_paths),...
        char(TCseg_paths), char(TPM_paths), param.fwhm_tws,param.l_TC, pth_out);
end

%% Apply tissue-specific smoothing compensated (TSPOON)
if flag.tspoon
    disp('Executing tissue-specific smoothing compensated (TSPOON)...');
    
    % Automatically create a new smoothing-specific derivatives directory for the subject
    pth_out = regexprep(spm_file(MPM_paths(1,:), 'fpath'), param.outDerivName, 'AJ-TSPOON');
    if ~exist(pth_out, 'dir')
        mkdir(pth_out);
    else
        % If the directory already exists, delete all the files inside it
        file_list = dir(fullfile(pth_out, '*'));
        for k = 1:length(file_list)
            % Ignore the special directories '.' and '..'
            if ~file_list(k).isdir
                delete(fullfile(pth_out, file_list(k).name));
            end
        end
    end
    
    [tspoon_paths, sMask_paths] = aj_proc_MPMTPSOON(char(MPM_paths), ...
        char(TCseg_paths), param.fwhm_tspoon,param.l_TC, pth_out);
end

end
