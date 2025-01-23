function [ggsP_GmWmCsf, ttwsP_signal, ttosP_signal] = aj_reorganize_signals_by_tissue(gsP_signal, twsP_signal, tosP_signal, dim, nb_tissue)
% Toy function to reorganize signals per tissue type rather than per subject
% 
% INPUT
% gsP_signal: Matrix of smoothed signals per subject (nb_tissue x ...).
% twsP_signal: Matrix of tissue-weighted signals (nb_tissue x ...).
% tosP_signal: Matrix of other signals (nb_tissue x ...).
% dim: Dimension of the data (1 for 1D, 2 for 2D, 3 for 3D)
% nb_tissue: Number of tissue types.
%
% OUTPUT
% ggsP_GmWmCsf: Cell array organized by tissue type for gsP_signal.
% ttwsP_signal: Cell array organized by tissue type for twsP_signal.
% ttosP_signal: Cell array organized by tissue type for tosP_signal.
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------
% Initialize cell arrays for storing results organized by tissue type
ggsP_GmWmCsf = cell(nb_tissue, 1);
ttwsP_signal = cell(nb_tissue, 1);
ttosP_signal = cell(nb_tissue, 1);

% Check the dimension and reorganize accordingly
switch dim
    case 1  % 1D case
        for jj = 1:nb_tissue
            ggsP_GmWmCsf{jj} = gsP_signal(jj, :);  % Store signal for tissue type jj
            ttwsP_signal{jj} = twsP_signal(jj, :);  % Store twsP_signal for tissue type jj
            ttosP_signal{jj} = tosP_signal(jj, :);  % Store tosP_signal for tissue type jj
        end

    case 2  % 2D case
        for jj = 1:nb_tissue
            ggsP_GmWmCsf{jj} = gsP_signal(jj, :, :);  % Store signal for tissue type jj in 2D
            ttwsP_signal{jj} = twsP_signal(jj, :, :);  % Store twsP_signal for tissue type jj in 2D
            ttosP_signal{jj} = tosP_signal(jj, :, :);  % Store tosP_signal for tissue type jj in 2D
        end

    case 3  % 3D case
        for jj = 1:nb_tissue
            ggsP_GmWmCsf{jj} = gsP_signal(jj, :, :, :);  % Store signal for tissue type jj in 3D
            ttwsP_signal{jj} = twsP_signal(jj, :, :, :);  % Store twsP_signal for tissue type jj in 3D
            ttosP_signal{jj} = tosP_signal(jj, :, :, :);  % Store tosP_signal for tissue type jj in 3D
        end

    otherwise
        error('Invalid dimension. dim must be 1, 2, or 3.');
end
end
