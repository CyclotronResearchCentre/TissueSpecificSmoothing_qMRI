function selected_files = aj_select_smoofiles(tmp_unzip_dir, smoo, smoo_dir, smoo_approach, qmetric, TC)

% Select files for each qMRI metric and each TC for all subjects
if smoo==1
    selected_files = cellstr(spm_select('FPListRec',smoo_dir,sprintf('^.*%ssmo_%s.*\\.nii$',TC,qmetric)));
elseif smoo==4
    pattern = sprintf('^susan_Mask_%s.*%s.*\\.nii(\\.gz)?$', TC, qmetric);
    file_list = cellstr(spm_select('FPListRec', smoo_dir, pattern));

    % Unzip files if needed
    unzipped_files = cell(size(file_list));

    for f = 1:numel(file_list)
        [~,base,ext] = fileparts(file_list{f});

        if strcmp(ext,'.gz')
            % full filename was .nii.gz → remove .gz
            nii_out = fullfile(tmp_unzip_dir, base);
            if ~exist(nii_out,'file')
                fprintf('Unzipping %s\n', file_list{f});
                gunzip(file_list{f}, tmp_unzip_dir);
            end
            unzipped_files{f} = nii_out;

        else
            % already .nii
            unzipped_files{f} = file_list{f};
        end
    end

    selected_files = unzipped_files;

else
    selected_files = cellstr(spm_select('FPListRec',smoo_dir,sprintf('^%s_%s.*%s.*\\.nii$',smoo_approach,TC,qmetric))); 
end

end