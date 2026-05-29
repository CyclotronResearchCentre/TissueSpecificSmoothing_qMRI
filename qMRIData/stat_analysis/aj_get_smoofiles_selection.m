function selected_files = aj_get_smoofiles_selection(smoo, smoo_dir, smoo_approaches, qmetric, TC)
% Select files for each qMRI metric and each TC for all subjects
if smoo==1
    pattern = sprintf('^.*%ssmo_%s.*\\.nii$',TC,qmetric);
elseif smoo==4
    pattern = sprintf('^susan_Mask_%s.*%s.*\\.nii(\\.gz)?$', TC, qmetric);
else
    pattern = sprintf('^%s_%s.*%s.*\\.nii$',smoo_approaches{smoo},TC,qmetric);
end

selected_files = cellstr(spm_select('FPListRec',smoo_dir,pattern)); 

end