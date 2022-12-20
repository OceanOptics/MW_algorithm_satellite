# MW_algorithm_satellite
Matlab code for the retrievals of SPM and its associated uncertainties adjusted to be applied to any satellite sensor

This module provides a package of Matlab functions to derived SPM and associated uncertianties fromsatellite remote sensing reflectance at any spectral resolution or satellite sensor. This is likely not the fastest implementation possible for the MW algorithm. However, it is easy to use and handles the processing of either multispectral and hyperspectral satellite remote sensing reflectance data. 

# Setting the input data 

The function 'RW2SPM.m' calls all the input information necessary and the sequence of functions to run the scripts.

Because radiometric satellite measurements are most likely to be made available as water-leaving reflectance (RW) we suggest the use of the function 'RW2rrs.m' to convert the measured water-leaving reflectance (RW) to below water remote sensing reflectance (rrs). 

In sequence the inputs necessary for for the main function 'MW_algorithm_sat.m' are the satellite band specific signal-to-noise (SNR) and bandwidhts, also the ranges of shape parameters (S and Y) and the mass-specific coefficientes (a_nap at 443 nm, a_nap at 750 nm, and b_bp at 700 nm), the water temperature corresponding to when RW was measured, and the saturation threshold og choice (Q_filter)

# Questions and Suggestions

For questions regarding the script implementation or to suggest changes to improve its functionality, please contact Juliana Tavora at juliana.tavora@maine.edu; j.tavora@utwente.nl
