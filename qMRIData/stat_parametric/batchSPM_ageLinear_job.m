%-----------------------------------------------------------------------
% Job saved on 20-Nov-2024 18:04:28 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
%% Change spm default values
spm_get_defaults('stats.fmri.ufp',0.5); % change nothing for qMRI except the computational speed
spm_get_defaults('stats.rft.nonstat',1); % set to 1 the Random Field Theory if assuming non stationary smoothing

disp(['Current stats.fmri.ufp: ', num2str(spm_get_defaults('stats.fmri.ufp'))]);
disp(['Current stats.rft.nonstat: ', num2str(spm_get_defaults('stats.rft.nonstat'))]);

%% GLM from Callaghan et al. (2014)
matlabbatch{1}.spm.stats.factorial_design.dir = '<UNDEFINED>';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.scans = '<UNDEFINED>';
load('E:\Master_Thesis\Data\BIDS_AgingData\reg_age.mat');
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).c = age;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).cname = 'Age';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(1).iCC = 1;
load('E:\Master_Thesis\Data\BIDS_AgingData\reg_sex.mat');
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(2).c = sex;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(2).cname = 'Sex';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(2).iCC = 1;
load('E:\Master_Thesis\Data\BIDS_AgingData\reg_TIV.mat');
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(3).c = TIV;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(3).cname = 'TIV';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(3).iCC = 1;
load('E:\Master_Thesis\Data\BIDS_AgingData\reg_scanner.mat');
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(4).c = scanner;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(4).cname = 'Scanner';
matlabbatch{1}.spm.stats.factorial_design.des.mreg.mcov(4).iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.des.mreg.incint = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = '<UNDEFINED>';
matlabbatch{1}.spm.stats.factorial_design.globalc.g_mean = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Increase with age';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [0 1 0 0 0]; % intercept in the first column
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Decrease with age';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 -1 0 0 0];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{3}.fcon.name = 'Depends on age';
matlabbatch{3}.spm.stats.con.consess{3}.fcon.weights = [0 1 0 0 0];
matlabbatch{3}.spm.stats.con.consess{3}.fcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;
