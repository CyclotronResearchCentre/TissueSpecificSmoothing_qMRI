% Function to create the regressors files if needed

function paths_reg = aj_compute_regfile(ds_dir)
tsv_path = fullfile(ds_dir, 'participants.tsv');

% Define the output file paths
reg_age_path = fullfile(ds_dir, 'reg_age.mat');
reg_sex_path = fullfile(ds_dir, 'reg_sex.mat');
reg_TIV_path = fullfile(ds_dir, 'reg_TIV.mat');
reg_scanner_path = fullfile(ds_dir, 'reg_scanner.mat');

% Check if all regressor files already exist
if exist(reg_age_path, 'file') && exist(reg_sex_path, 'file') && exist(reg_TIV_path, 'file') && exist(reg_scanner_path, 'file')
    fprintf('All regressor files already exist. Skipping creation.\n');
else
    fprintf('One or more regressor files do not exist. Creating regressor files...\n');

    % Read the .tsv file
    opts = detectImportOptions(tsv_path, 'FileType', 'text');
    data = readtable(tsv_path, opts);

    % Convert categorical variables (sex, scanner) to boolean
    sex = double(strcmp(data.sex, 'M')); % Male=1, Female=0
    scanner = double(strcmp(data.scanner, 'trio')); % Trio=1, Quatro=0
    age = data.age;
    TIV = data.TIV; % Total Intracranial Volume

    % Save the regressor files
    save(reg_age_path, 'age');
    save(reg_sex_path, 'sex');
    save(reg_TIV_path, 'TIV');
    save(reg_scanner_path, 'scanner');

    fprintf('Regressor files created successfully.\n');
end

paths_reg.age = reg_age_path;
paths_reg.sex = reg_sex_path;
paths_reg.TIV = reg_TIV_path;
paths_reg.scanner = reg_scanner_path;

end