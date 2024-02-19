% NCP in moles of oxygen per square metre are taken from the biology
% contribution term calculated in script 11, converting to per cubic metre
% by dividing at each timestep by the MLD.

mooring_data.ncp_O2_molm2hr = mooring_data.dox2_biology_molm2;

mooring_data.ncp_O2_molm3hr = mooring_data.dox2_biology_molm2./mooring_data.mld_smooth;

mooring_data.ncp_O2_umolkghr = 1E6*mooring_data.ncp_O2_molm3hr./mooring_data.density_kgm3;
 
% This is converted into milligrams of carbon per m^2 per hour, dividing
% by 1E6 and a default Redfield ratio of 1.45 (Anderson and Sarmiento, 1994), 
% and multiplying by the molar mass of Carbon, and 1000.
 
mooring_data.ncp_C_mgm2hr = mooring_data.ncp_O2_molm2hr/constants.ncp_oxygen2carbon*constants.atomic_mass_C*1000;


% The NCP estimates are cummulatively summed to arrive at cummulative time
% series, with the total NCP for the timeseries being the final value of
% each vector

mooring_data.ncp_C_mgm2_cumsum = cumsum(mooring_data.ncp_C_mgm2hr,'omitnan');

mooring_data.ncp_O2_molm2_cumsum = cumsum(mooring_data.ncp_O2_molm2hr,'omitnan');
