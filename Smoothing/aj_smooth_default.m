function [param, flag] = aj_smooth_default()
%--------------------------------------------------------------------------
% Function to set default parameters and flags before computing smoothing
%--------------------------------------------------------------------------
% Copyright (C) 2017 Cyclotron Research Centre
% Written by A.J.
% Cyclotron Research Centre, University of Liege, Belgium
%--------------------------------------------------------------------------

% Parameters Definitions
param.fwhm_gs = 5;          % Kernel width for Gaussian Smoothing (GS)
param.fwhm_tws = 5;         % Kernel width for Tissue-Weighted Smoothing (TWS)
param.fwhm_tspoon = 5;      % Kernel width for Tissue-SPecific smOOthing compeNsated (TSPOON)

param.l_TC = (1:2);

% Flag: Used to manage additional behaviors
flag.plot_fig = false;           % Flag to plot figures (set to true for plotting)

flag.gaussian = false;           % Flag to compute gaussian smoothing
flag.tws = false;                % Flag to compute Tissue-Weighted Smoothing
flag.tspoon = true;             % Flag to compute Tissue-SPecific smOOthing compeNsated

end