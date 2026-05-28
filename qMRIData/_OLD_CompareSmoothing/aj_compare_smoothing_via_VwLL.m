% aj_compare_smoothing_via_VwLL.m
% ---------------------------------------------------------
% Robust comparison of smoothing approaches using voxelwise GLM likelihood
% Author: Antoine Jacquemin
% ---------------------------------------------------------
%% ================= SETTINGS =================
clear; clc;

addpath('C:\Users\antoi\Documents\JOBS\PhD\MaterThesis\scripts\spm12');
spm('Defaults','fMRI'); spm_jobman('initcfg');

combination_names = {
    'MTsat_GM'
    'MTsat_WM'
    'PDmap_GM'
    'PDmap_WM'
    'R1map_GM'
    'R1map_WM'
    'R2starmap_GM'
    'R2starmap_WM'
};

% combination_names = {
%     'MTsat_GM'
% };

basepaths = {
    'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-SUSAN_GLM_rft1\'
    'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TSPOON_GLM_rft1\'
    'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\AJ-TWS_GLM_rft1\'
};

out_root = 'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\LL_AIC_BIC_results\parametric_rft1';

beta_names = {'beta_0001.nii','beta_0002.nii','beta_0003.nii','beta_0004.nii','beta_0005.nii'};
ResMS_name = 'ResMS.nii';

% Correction of scans' paths
old_root = 'C:\Users\lucad\Documents\smoothing\data\qMRI_AgingCallaghan\derivatives\';
new_root = 'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\';

chunk_vox = 60000;
sigma2_floor = 1e-8;
verbose = true;

%% ================= MAIN LOOP =================
for c = 1:numel(combination_names)

    comb = combination_names{c};
    fprintf('\n===== %s =====\n', comb);

    %% ---------- MASK ----------
    if contains(comb,'WM')
        mask_fn = 'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\mask_WTA_WM_mask.nii';
    else
        mask_fn = 'E:\Master_Thesis\Data\BIDS_AgingData\derivatives\mask_WTA_GM_mask.nii';
    end

    Vmask = spm_vol(mask_fn);
    mask = spm_read_vols(Vmask) > 0;
    vox_idx = find(mask);
    nvox = numel(vox_idx);

    %% ---------- LOAD METHODS ----------
    M = numel(basepaths);
    X_all = cell(M,1);
    scans = cell(M,1);
    betas = cell(M,1);
    ResMS = cell(M,1);
    df_all = zeros(M,1);

    for m = 1:M
        spm_path = fullfile(basepaths{m}, comb, 'SPM.mat');
        S = load(spm_path); SPM = S.SPM;

        X_all{m} = SPM.xX.X;
        df_all(m) = SPM.xX.erdf;
        scans{m} = SPM.xY.P;

        % beta
        for k = 1:numel(beta_names)
            betas{m}{k} = spm_vol(fullfile(basepaths{m}, comb, beta_names{k}));
        end

        ResMS{m} = spm_vol(fullfile(basepaths{m}, comb, ResMS_name));
    end
    
    % Update scan paths
    scans = correct_scan_path(old_root,new_root,scans);

    Nsub = size(X_all{1},1);

    %% ---------- OUTPUT ----------
    outdir = fullfile(out_root, comb);
    if ~exist(outdir,'dir'), mkdir(outdir); end

    Vtempl = spm_vol(scans{1}{1});

    % maps
    LL_maps = init_maps(Vtempl, outdir, M, 'LL');
    best_map = init_single_map(Vtempl, outdir, 'best_LL');

    % NEW
    deltaLL_map = init_maps(Vtempl, outdir, M, 'deltaLL');
    softmax_map = init_maps(Vtempl, outdir, M, 'softmax');

    %% ---------- CHUNKING ----------
    blocks = chunk_indices(vox_idx, chunk_vox);

    %% ================= CORE =================
    for b = 1:numel(blocks)

        vox = blocks{b};
        nb = numel(vox);

        LL = zeros(nb,M);

        for m = 1:M

            X = X_all{m};
            p = size(X,2);

            % ---- betas ----
            B = zeros(nb,p);
            for k = 1:p
                vol = betas{m}{k};
                img = spm_read_vols(vol);          % lire tout le volume
                B(:,k) = img(vox);                 % indexer avec vox
            end

            % ---- data ----
            Y = zeros(nb,Nsub);
            for s = 1:Nsub
                vol = spm_vol(scans{m}{s});
                img = spm_read_vols(vol);          % lire tout le volume
                Y(:,s) = img(vox);                 % indexer avec vox
            end

            % ---- prediction ----
            mu = (X * B')';

            % ---- variance (CRUCIAL FIX) ----
            sigma2 = spm_get_data(ResMS{m}, vox);
            sigma2 = sigma2 * (df_all(m) / Nsub);
            sigma2(sigma2 < sigma2_floor) = sigma2_floor;

            % ---- LL ----
            resid2 = (Y - mu).^2;
            LL(:,m) = -0.5 * ( ...
                Nsub .* log(2*pi*sigma2) + ...
                sum(resid2,2) ./ sigma2 ...
            );
        end

        %% ====== METRICS ======

        % 1. BEST LL
        [~, best_idx] = max(LL,[],2);

        % 2. DELTA LL (relative evidence)
        maxLL = max(LL,[],2);
        deltaLL = LL - maxLL;

        % 3. SOFTMAX (posterior prob)
        expLL = exp(deltaLL); % stable
        softmax = expLL ./ sum(expLL,2);

        %% ====== WRITE ======
        write_maps(LL_maps, LL, vox);
        write_maps(deltaLL_map, deltaLL, vox);
        write_maps(softmax_map, softmax, vox);
        write_single_map(best_map, best_idx, vox);

        fprintf('Block %d/%d done\n', b, numel(blocks));
    end

    fprintf('DONE %s\n', comb);
end

%% ================= FUNCTIONS =================

function maps = init_maps(Vtempl, outdir, M, prefix)
    maps = cell(M,1);
    for m = 1:M
        V = Vtempl;
        V.fname = fullfile(outdir, sprintf('%s_m%d.nii',prefix,m));
        spm_write_vol(V, zeros(V.dim));
        maps{m} = V.fname;
    end
end

function V = init_single_map(Vtempl, outdir, name)
    V = Vtempl;
    V.fname = fullfile(outdir, [name '.nii']);
    spm_write_vol(V, zeros(V.dim));
end

function blocks = chunk_indices(idx, chunk_size)
    n = numel(idx);
    nb = ceil(n/chunk_size);
    blocks = cell(nb,1);
    for i = 1:nb
        blocks{i} = idx((i-1)*chunk_size+1 : min(i*chunk_size,n));
    end
end

function write_maps(map_files, data, vox)
    for m = 1:numel(map_files)
        V = spm_vol(map_files{m});
        img = spm_read_vols(V);
        img(vox) = data(:,m);
        spm_write_vol(V,img);
    end
end

function write_single_map(Vfile, data, vox)
    V = spm_vol(Vfile.fname);
    img = spm_read_vols(V);
    img(vox) = data;
    spm_write_vol(V,img);
end

function scans = correct_scan_path(old_root,new_root,scans)

    for m = 1:numel(scans)

        % --- convertir en cell array si nécessaire ---
        if ischar(scans{m})
            scans{m} = cellstr(scans{m});
        elseif isnumeric(scans{m})
            error('scans{%d} is numeric → unexpected format', m);
        end

        % --- remplacer les paths ---
        scans{m} = cellfun(@(p) strrep(p, old_root, new_root), ...
                           scans{m}, 'UniformOutput', false);
    end
end