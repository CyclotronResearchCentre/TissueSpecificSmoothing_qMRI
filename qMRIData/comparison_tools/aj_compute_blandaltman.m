function aj_compute_blandaltman(LL_dir, mask_fn, methodA, methodB, flag)
% === Load LL maps ===
LL_A_fn = fullfile(LL_dir, sprintf('LL_method%d.nii', methodA));
LL_B_fn = fullfile(LL_dir, sprintf('LL_method%d.nii', methodB));

V_A = spm_vol(LL_A_fn);
V_B = spm_vol(LL_B_fn);

LL_A = spm_read_vols(V_A);
LL_B = spm_read_vols(V_B);

% Output folder for Bland–Altman
outdir = fullfile(LL_dir,'BlandAltman_LL');
if ~exist(outdir,'dir'), mkdir(outdir); end

% === Load tissue mask (GM or WM) ===
Vmask = spm_vol(mask_fn);
mask = spm_read_vols(Vmask) > 0;

% Extract only valid voxels
LL_A_vals = LL_A(mask);
LL_B_vals = LL_B(mask);


% === Remove NaNs or inf (security) ===
valid = ~(isnan(LL_A_vals) | isnan(LL_B_vals) | isinf(LL_A_vals) | isinf(LL_B_vals));
LL_A_vals = LL_A_vals(valid);
LL_B_vals = LL_B_vals(valid);


% === Bland–Altman (using your aj_BlandAltman function) ===

% Define method names
method_names = {'SUSANs', 'gTSPOON', 'TWS'};

% Create plot title with method names
plot_title = sprintf('Bland–Altman LL: %s vs %s', method_names{methodA}, method_names{methodB});

% Create output filename
out_eps = fullfile(outdir, sprintf('LL_%s_%s.png', method_names{methodA}, method_names{methodB}));

% Add "_fit" to the output filename if fitting is enabled
if isfield(flag, 'fitting') && flag.fitting == 1
    [filepath, name, ext] = fileparts(out_eps);
    out_eps = fullfile(filepath, [name, '_fit', ext]);
end

% plot_title = sprintf('Bland–Altman LL: method %d vs %d', methodA, methodB);
% out_eps    = fullfile(outdir, sprintf('LL_method%d_method%d.png', methodA, methodB));

[mean_diff, std_diff] = aj_BlandAltman(LL_A_vals, LL_B_vals, flag, out_eps, plot_title);
% Note: Mean difference (meth1 - meth2) in aj_BlandAltman

% === Print summary ===
fprintf('\n==== Bland–Altman summary (LL method %d vs %d) ====\n', methodA, methodB);
fprintf('Mean difference (A–B): %.4f\n', mean_diff);
fprintf('Std of difference     : %.4f\n', std_diff);
fprintf('Limits of agreement   : [%.4f , %.4f]\n', ...
         mean_diff - 1.96*std_diff, mean_diff + 1.96*std_diff);
fprintf('Plot saved at: %s\n', out_eps);
end