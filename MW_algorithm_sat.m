function [SPM_mw, err_mw] = MW_algorithm(nm, std, rrs, U, ...
    IOP_matrix, a_sw, Q_filter, SNR, dim,temp)
% MW algorithm and for SPM and associated uncertainty estimates
%
% Juliana Tavora, University of Maine, 2020
%
% See the following publication for details on the method:
% Tavora, J, et al., An algorithm to estimate Suspended Particulate Matter
% concentrations and associated uncertainties from Remote Sensing Reflectance
% in Coastal Environments
%
% INPUTS:
%
% nm         -  wavelengths associated with measured Remote sensing reflectance
% std        -  standard deviation of measared below water rrs(nm) 
% rrs        -  measured below water rrs(nm)
% U          -  function of IOPs by Gordon et al (1988) 
% L          -  constants to calculate U by Gordon et al (1988)
% S          -  exponent of exponential spectral shape for a_NAP :[0.006:0.001:0.014]
% Y          -  exponent of power-law spectral shape for bbp : [0:0.15:1.8]
% a_nap443   -  mass specific absorption at the reference at 443 nm : [0.01:0.01:0.06]
% a_nap750   -  mass specific absorption at the reference at 750 nm :[0.013:0.001:0.015]
% bbp700     -  mass specific backscattering at 700 nm : [0.002:0.001:0.021]
% temp       -  temperature of the water in Celsius which remote sensing
%               reflectance was measured
% Q_filter   -  saturation threshold 
%
% This function calls other 3 functions:  
%   1) 'asw_correction' Water absorption for temperature correction 
%   2) 'weight_assess' for estimates of maximum uncertainties and
%       establishing the weight used for final SPM and uncertainty estimates
%   3) 'IOP_sel' constrains combinations of IOPs and returns a table with
%   the 20 most common IOPs for each measurement
%
% OUTPUTS:
%
% SPM_mw    -  SPM estimates
% err_mw    -  uncertainties associated with the method in gm^-3         
% IOP_table -  table of the 20 most common combinations of IOPs. The
%              following code line can be added to estimate median and ranges
%              of IOP
%
% The SPM retrivals and uncertainties can be applied to ocean color 
% remote sensing data such as Landsat8/OLI, MODIS/aqua and others. 
% Uncertainties calculated here do not take into account the uncertainties 
% regarding atmospheric correction uncertainties (see section uncertainties,
% reference above).

%-------------------------------------------------------------------------%
L = [0.0949,0.0794];
%--------------------------------------------------------------------------%
%% SPM calculation

Q = NaN(length(nm),(dim),length(rrs(:,1)));


U_model = (squeeze(IOP_matrix(2,:,:)) ./ (squeeze(IOP_matrix(1,:,:)) + squeeze(IOP_matrix(2,:,:))));


ii=1;
for i=1:length(rrs(:,1))
    %----------------------------------------------------------------------------------------------------------------------------------------------------%
    % SPM estimates
    
    if length(temp) ==1
        SPM_model(:,:,ii) = (1./squeeze(IOP_matrix(2,:,:)) .* (a_sw ./  ( (((1-U(i,:))./U(i,:))' - (squeeze(IOP_matrix(1,:,:))./squeeze(IOP_matrix(2,:,:)))))));
    else
        SPM_model(:,:,ii) = (1./squeeze(IOP_matrix(2,:,:)) .* (squeeze(a_sw(:,:,i)) ./  ( (((1-U(i,:))./U(i,:))' - (squeeze(IOP_matrix(1,:,:))./squeeze(IOP_matrix(2,:,:)))))));
    end
    %----------------------------------------------------------------------------------------------------------------------------------------------------%
    % saturation estimates
    Q(:,:,ii) = U(i,:)' ./ U_model;
    
    %----------------------------------------------------------------------------------------------------------------------------------------------------%
    fprintf('%0.0f/%0.0f\n',i,length(rrs(:,1)));
    
    ii=ii+1;
    
end

SPM_model(SPM_model < 0 | Q < 0 | Q > Q_filter) = NaN;

%---------------------------------------------------------------------------------------------------------------------------------------------------------%
%% Solutions

% ----------------------------------------------------------------------- %
% Determing Degrees of Freedom of rrs spectra (i.e., number of satellite
% bands used)
M = length(nm);

% ----------------------------------------------------------------------- %
% Weighting estimates 

[Weight] = weight_asses(std, rrs, SPM_model, IOP_matrix, L, U, nm, Q, Q_filter, SNR);
 W = repmat(Weight,1,1,(dim));
 W = permute(W, [1 3 2]);

% filtering results for saturation 
W(SPM_model < 0 | Q < 0 | Q > Q_filter)=NaN;
Q(SPM_model < 0 | Q < 0 | Q > Q_filter)=NaN;
coefs_matrix(:, SPM_model < 0 | Q < 0 | Q > Q_filter) = NaN;
shape_matrix(:, SPM_model < 0 | Q < 0 | Q > Q_filter) = NaN;

% ----------------------------------------------------------------------- %
% SPM retrieval

SPM_mw = (squeeze(nansum(nansum(SPM_model.*W),2)) ./ squeeze(nansum(nansum(W,2))))';

SPM_84 = squeeze(prctile(squeeze(SPM_model),84,2));
SPM_16 = squeeze(prctile(squeeze(SPM_model),16,2));

SPM_max = ((nansum(SPM_84.*Weight)) ./ (nansum(Weight))).*(1./sqrt(M));
SPM_min = ((nansum(SPM_16.*Weight)) ./ (nansum(Weight))).*(1./sqrt(M));

err_mw = ((SPM_max - SPM_min)./2);


end
