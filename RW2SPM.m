
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
[U, rrs, nm]           = RW2rrs(RW, nm, L);
[IOP_matrix,a_sw,dim]  = gen_IOP_array(rrs, U, std, bandwidth, SNR, nm, temp, S, Y, ap443, ap750, bbp700);
[SPM_mw, err_mw]       = MW_algorithm_sat(nm, std, rrs, U, IOP_mtrix, asw, Q_filter, SNR, dim, temp);
        
        
%% OUTPUTS:

%SPM_modeled:
SPM_mw;

%SPM_uncertainty (in gm^-3):
err_mw;




