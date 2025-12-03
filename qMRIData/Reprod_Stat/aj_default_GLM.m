function [paths] = aj_default_GLM(flags)
%% Define a default values for flags
if ~isfield(flags, 'users'), flags.users = -1; end

%% Do the job
% Paths to access to scripts, data and spm toolbox
if flags.users == 0 % Personal PC
    paths.script_dir = 'C:\Users\antoi\Documents\JOBS\PhD\MaterThesis\scripts\TissueSpecificSmoothing\qMRIData\Reprod_Stat';
    paths.ds_dir = 'E:\Master_Thesis\Data\BIDS_AgingData';
    addpath("C:\Users\antoi\Documents\JOBS\PhD\MaterThesis\scripts\spm12");
elseif flags.users == 1 % Desktop at CRC
    paths.script_dir = 'C:\Users\lucad\Documents\smoothing\repo\TissueSpecificSmoothing\qMRIData\Reprod_Stat';
    paths.ds_dir = 'C:\Users\lucad\Documents\smoothing\data\qMRI_AgingCallaghan';
    addpath("C:\Users\lucad\Documents\smoothing\repo\spm12");
elseif flags.users == 99 % Skip setting paths
    % Do nothing, just continue without modifying param.script_dir and param.ds_dir
else % non existing user
    warning('Unexpected user flag: %d. Using default paths or none.', flags.users);
    paths.script_dir = '';
    paths.ds_dir = '';
end
end