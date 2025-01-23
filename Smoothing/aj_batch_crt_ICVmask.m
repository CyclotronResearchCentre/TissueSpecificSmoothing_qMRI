function outpath_ICVmask = aj_batch_crt_ICVmask()
%--------------------------------------------------------------------------
% Function to create ICV mask amoung the 16 subjects of ds000117. 
%
% OUTPUT
% outpath_ICVmask:  full path to the ICVmask that represents areas where 
%                   the probability of being part of the intracranial 
%                   volume (ICV) is greater than 50%.
%
% BATCH PROCESS
% 1. Mean mwc1, mwc2 and mwc3 amoung all subjects using the formula
%    sum(X)/(size(X,1)/3) -> mean_mwc123.nii
% 2. Smooth mean_mwc123.nii with fwhm = [2 2 2] -> smean_mwc123.nii
% 3. Apply a threshold of 0.5 using spm_imcalc -> ICVmask_pop.nii
%--------------------------------------------------------------------------
% FUTURE : manage dir dependencies (here + in batch)
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------

nsub = 16;
dir = 'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\derivatives\preprocessing';

inputs = cell(nsub,1);
for s = 1:nsub
    inputs{s} = cellstr(spm_select('ExtFPList', ...
        fullfile(dir, sprintf('sub-%02d', s), 'ses-mri', 'anat'), ...
        '^mwc.*\.nii$', Inf));
end

jobfile = {fullfile(pwd, 'aj_batch_crt_ICVmask_job.m')};

spm('defaults', 'FMRI');
spm_jobman('run', jobfile, inputs{:});

% Delete unuseful files
% delete(fullfile(dir,'mean_mwc123.nii'));
% delete(fullfile(dir,'smean_mwc123.nii'));

% Full path to output
outpath_ICVmask = fullfile(dir,'ICVmask_pop.nii');
% spm_image('Display',outpath_ICVmask);

% a=spm_read_vols(spm_vol(fullfile(dir,'smean_mwc123.nii')));
% b=spm_read_vols(spm_vol(fullfile(dir,'ICVmask_pop.nii')));

end
