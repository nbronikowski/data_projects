% Calculate required Air-Sea N2 exchange to match N2 record

%%
% We convert the gas transfer velocities into metres per second.
mooring_data.Gc_O2_ms = mooring_data.Gc_O2_cmhr./(100*3600);

mooring_data.Gc_N2_ms = mooring_data.Gc_N2_cmhr./(100*3600);


 
% We use Gc to calculate the rate of exchange of gas between atmosphere and
% water in moles per m^2 per s, as described in Emerson et al. 2008, eqns 4
% and 5.

mooring_data.GE_N2 = mooring_data.Gc_N2_ms.*(mooring_data.N2_molm3 - (mooring_data.atmosphericpress_Pa/constants.atm_in_Pa).*mooring_data.N2sol_molm3);
  
% In order to estiamte gas injections due to bubbles, we need to estimate
% molecular diffusion.

%    O2 and N2 freshwater values from  Ferrell and Himmelblau, 1967.
%       "Diffusion coefficients of nitrogen and oxygen in water"
%       J. Chem. Eng. Data, 12(1), 111-115, doi: 10.1021/je60032a036.
%    Correction for salinity is based on Jahne's observed average 4.9% decrease in 
%       diffusivity for H2 and He in 35.5 ppt NaCl solution
% Roberta Hamme 2013

mooring_data.MDiff_O2 = (constants.Ac_O2*exp(constants.Ea_O2./(constants.gas_constant*1E-3*(mooring_data.temp_C+273.15))))*(1E-5).*(1 - constants.Jahne_sal_co*mooring_data.psal_PSU/35.5); % in cm^2 per s

mooring_data.MDiff_N2 = (constants.Ac_N2*exp(constants.Ea_N2./(constants.gas_constant*1E-3*(mooring_data.temp_C+273.15))))*(1E-5).*(1 - constants.Jahne_sal_co*mooring_data.psal_PSU/35.5); % in cm^2 per s


% We smooth the mixed layer depth time series, differentiate wrt time, and
% set times of decreasing depth to 0, for the purpose of calculating
% entrainment, as we are only interested when the MLDP deepens, bringing
% oxygen into the layer. 

mooring_data.mld_smooth = smooth(mooring_data.mld_m,mld_smooth_span);

mooring_data.mld_diff = diff(mooring_data.mld_smooth);

mooring_data.mld_increasing = find(mooring_data.mld_diff<=0);

mooring_data.mld_diff(mooring_data.mld_increasing) = 0;

% We preallocate zero arrays used in the for loop

[mooring_data.N2_gas_ex_molm3, mooring_data.N2_bubbles_molm3, mooring_data.Vinj] = deal(zeros(size(mooring_data.time)));

for i = 2:length(mooring_data.time)
    
    % This calculates the N2 gas exchange from the previous time stamp using
    % the calculate gas exchange rate from that time stamp, and the MLDP of
    % the current time stamp (as we are calculating the gas that will have
    % exchanged in the previous time stamp into the current time stamp).
    
    mooring_data.N2_gas_ex_molm3(i-1) = -3600*mooring_data.GE_N2(i-1)/mooring_data.mld_smooth(i);
    
    % This subtracts the calculated gas exchange from the measured change, to
    % determine the remainder, which is attribute to bubbles.
    
    mooring_data.N2_bubbles_molm3(i-1) = mooring_data.N2_molm3(i) - mooring_data.N2_molm3(i-1) - mooring_data.N2_gas_ex_molm3(i-1);
    
    % This rearranges eqn 9 from Emerson et al. 2008 to solve for Vinj, in units of moles per square metre per second
    mooring_data.Vinj(i-1) = mooring_data.N2_bubbles_molm3(i-1)/(3600*constants.mole_fraction_N2*(1+(1/bubble_beta)*((mooring_data.MDiff_N2(i-1)/mooring_data.MDiff_O2(i-1))^(0.5))*(mooring_data.N2sol_umolkg(i-1)/mooring_data.dox2_sol_umolkg(i-1)))/mooring_data.mld_smooth(i)); 
    
end

% As the calculation of at a timestamp Vinj requires the N2 value for the 
% timestamp ahead of it, we assign the final Vinj value as Nan.

mooring_data.Vinj(length(mooring_data.time)) = nan;







