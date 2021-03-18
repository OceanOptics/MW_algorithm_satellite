
%% script example for the MW algorithm

% MW algorithm and for SPM and associated uncertainty estimates
%
% Juliana Tavora, University of Maine, 2020
%
% See the following publication for details on the method:
% Tavora, J, et al., An algorithm to estimate Suspended Particulate Matter
% concentrations and associated uncertainties from Remote Sensing Reflectance
% in Coastal Environments
%

close all
clear all; 
clc

%-------------------------------------------------------------------------%
%% INPUTS:

L         = [0.0949,0.0794]; %Gordon et al., 1988 constants
Q_filter  = 0.5; % saturation threshold 

% shape parameters
S         = 0.006:0.001:0.014;
Y         = 0:0.15:1.8;

% mass-specific coefficients
a_nap443     = 0.01:0.01:0.06;
a_nap750     = 0.013:0.001:0.015; 
bbp700       = 0.002:0.001:0.021;

%-------------------------------------------------------------------------%
% Setting up satellite specifc inputs required on the MW_algorithm_sat.m

SNR       = []; % satellite sensor and band specific signal-to-noise 
bandwidth = []; % satellite sensor and band specific bandwidth
temp      = []; % water temperature relative to when RW data were acquired

RW  = []; % water reflectances
std = []; % standard deviation of below water remote sensing reflectances
nm  = []; % central wavelength of satellite bands 

%% Water-leaving reflectance to below water remote sensing reflectance
[U, rrs, nm] = RW2rrs(RW, nm, L);
[SPM_mw, err_mw, IOP_table,SPM_max,SPM_min] = MW_algorithm_sat(nm, std, rrs, U, L, S, Y, a_nap443, bbp700, a_nap750, temp, Q_filter, SNR, bandwidth);

%% OUTPUTS:

%SPM_modeled:
SPM_mw;

%SPM_uncertainty (in gm^-3):
err_mw;

% constrained shape paremeters and mass-specific coeffs: 
IOP_S       = IOP_table(:,:,6);
IOP_Y       = IOP_table(:,:,5);
IOP_ap443   = IOP_table(:,:,3);
IOP_ap750   = IOP_table(:,:,4);
IOP_bbp700  = IOP_table(:,:,2);




