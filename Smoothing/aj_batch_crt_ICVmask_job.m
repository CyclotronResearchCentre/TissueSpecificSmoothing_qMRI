%-----------------------------------------------------------------------
% Job saved on 14-Nov-2024 09:52:50 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
%%
matlabbatch{1}.spm.util.imcalc.input = {'<UNDEFINED>'
                                        '<UNDEFINED>'
                                        '<UNDEFINED>'
                                        '<UNDEFINED>'
                                        '<UNDEFINED>'
                                        '<UNDEFINED>'
                                        '<UNDEFINED>'
                                        '<UNDEFINED>'
                                        '<UNDEFINED>'
                                        '<UNDEFINED>'
                                        '<UNDEFINED>'
                                        '<UNDEFINED>'
                                        '<UNDEFINED>'
                                        '<UNDEFINED>'
                                        '<UNDEFINED>'
                                        '<UNDEFINED>'
                                       };
%%
matlabbatch{1}.spm.util.imcalc.output = 'mean_mwc123.nii';
matlabbatch{1}.spm.util.imcalc.outdir = {'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\derivatives\preprocessing'};
matlabbatch{1}.spm.util.imcalc.expression = 'sum(X)/(size(X,1)/3)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Image Calculator: ImCalc Computed Image: mean_mwc123.nii', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{2}.spm.spatial.smooth.fwhm = [2 2 2];
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = 's';
matlabbatch{3}.spm.util.imcalc.input(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{3}.spm.util.imcalc.output = 'ICVmask_pop.nii';
matlabbatch{3}.spm.util.imcalc.outdir = {'C:\Users\antoi\Documents\master_thesis\MATLAB\ds000117\work_copies\openneuro.org\ds000117\derivatives\preprocessing'};
matlabbatch{3}.spm.util.imcalc.expression = 'i1>.5';
matlabbatch{3}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{3}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{3}.spm.util.imcalc.options.mask = 0;
matlabbatch{3}.spm.util.imcalc.options.interp = 1;
matlabbatch{3}.spm.util.imcalc.options.dtype = 2;
