% Function to rename GTSPOON files if needed

function aj_rename_gtspoon(smoo_dir)
% Check if there are GM files that need renaming
GM_files = spm_select('FPListRec', smoo_dir, sprintf('^tspoon1_sub*'));
need_renaming_GM = ~isempty(deblank(GM_files));

% Check if there are WM files that need renaming
WM_files = spm_select('FPListRec', smoo_dir, sprintf('^tspoon2_sub*'));
need_renaming_WM = ~isempty(deblank(WM_files));

% If any files need renaming, proceed with the renaming script
if need_renaming_GM || need_renaming_WM
    fprintf('Files found that need renaming. Starting renaming process...\n');

    % Rename GM files
    if need_renaming_GM
        nfiles = size(GM_files, 1);
        for i = 1:nfiles
            original_file = strtrim(GM_files(i, :));
            original_basename = spm_file(original_file, 'basename');
            new_basename = regexprep(original_basename, 'tspoon1_', 'TSPOON_GM_');
            new_file = fullfile(fileparts(original_file), [new_basename, '.nii']);

            movefile(original_file, new_file);
            fprintf('Renamed: %s -> %s\n', original_file, new_file);
        end
    end

    % Rename WM files
    if need_renaming_WM
        nfiles = size(WM_files, 1);
        for i = 1:nfiles
            original_file = strtrim(WM_files(i, :));
            original_basename = spm_file(original_file, 'basename');
            new_basename = regexprep(original_basename, 'tspoon2_', 'TSPOON_WM_');
            new_file = fullfile(fileparts(original_file), [new_basename, '.nii']);

            movefile(original_file, new_file);
            fprintf('Renamed: %s -> %s\n', original_file, new_file);
        end
    end
else
    fprintf('No files found that need renaming. Skipping renaming process.\n');
end
end

%% OLD script
% param.smoothName = 'AJ-TSPOON';
% ds_dir = 'D:\Master_Thesis\Data\BIDS_AgingData';
% smoo_dir = fullfile(ds_dir, 'derivatives', param.smoothName);
% 
% % Get list of GM files
% GM_files = spm_select('FPListRec', smoo_dir, sprintf('^tspoon1_sub*'));
% nfiles = size(GM_files, 1);
% 
% % Rename each file
% for i = 1:nfiles
%     original_file = strtrim(GM_files(i, :));
%     original_basename = spm_file(original_file, 'basename');
%     new_basename = regexprep(original_basename, 'tspoon1_', 'TSPOON_GM_');
%     new_file = spm_file(original_file, 'basename', new_basename);
%     movefile(original_file, new_file);
%     
%     fprintf('Renamed: %s -> %s\n', original_file, new_file);
% end
% 
% % Get list of WM files
% WM_files = spm_select('FPListRec', smoo_dir, sprintf('^tspoon2_sub*'));
% nfiles = size(WM_files, 1);
% 
% % Rename each file
% for i = 1:nfiles
%     original_file = strtrim(WM_files(i, :));
%     original_basename = spm_file(original_file, 'basename');
%     new_basename = regexprep(original_basename, 'tspoon2_', 'TSPOON_WM_');
%     new_file = spm_file(original_file, 'basename', new_basename);
%     movefile(original_file, new_file);
%     
%     fprintf('Renamed: %s -> %s\n', original_file, new_file);
% end