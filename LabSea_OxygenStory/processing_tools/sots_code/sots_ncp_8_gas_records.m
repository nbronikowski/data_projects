 % Generates Henry's law constants, nitrogen, and water vapour time
% series

%%
% Here we calculate N2 solubility in um/kg at 1atm from Hamme and Emerson
% 2004.
 
% We calculate a scaled temperature
mooring_data.temp_scaled = log((298.15 - mooring_data.temp_C)./(273.15 + mooring_data.temp_C));
 
% Eqn (1) of Hamme and Emerson 2004
mooring_data.N2sol_umolkg = exp(constants.A0_n2 + constants.A1_n2*mooring_data.temp_scaled + constants.A2_n2*mooring_data.temp_scaled.^2 + constants.A3_n2*mooring_data.temp_scaled.^3 + mooring_data.psal_PSU.*(constants.B0_n2 + constants.B1_n2*mooring_data.temp_scaled + constants.B2_n2*mooring_data.temp_scaled.^2));

%%
% Calculation of vapour pressure of seawater from Wagner and Pruss (2002)
% in Dickson, Sabine et al. (2007). 
 
% Calculate temperature in Kelvin and modified temperature for Chebyshev polynomial
mooring_data.temp_K = mooring_data.temp_C+273.15;
mooring_data.temp_mod = 1-mooring_data.temp_K./647.096;

% Calculate value of Wagner polynomial
mooring_data.Wagner = constants.wagner_a1*mooring_data.temp_mod +constants.wagner_a2*mooring_data.temp_mod.^1.5 +constants.wagner_a3*mooring_data.temp_mod.^3 +constants.wagner_a4*mooring_data.temp_mod.^3.5 +constants.wagner_a5*mooring_data.temp_mod.^4 +constants.wagner_a6*mooring_data.temp_mod.^7.5;

% Calculate vapour pressure of pure water in kiloPascals
mooring_data.vapour_0sal_kPa = exp(mooring_data.Wagner * 647.096 ./ mooring_data.temp_K) .* 22.064 * 1000;
 
% Correct vapour pressure for salinity
mooring_data.molality = 31.998 * mooring_data.psal_PSU ./(1e3-1.005*mooring_data.psal_PSU);
mooring_data.osmotic_coef = constants.millero_c0 +constants.millero_c1*(0.5*mooring_data.molality) +constants.millero_c2*(0.5*mooring_data.molality).^2 +constants.millero_c3*(0.5*mooring_data.molality).^3 +constants.millero_c4*(0.5*mooring_data.molality).^4;
mooring_data.vapour_press_kPa = mooring_data.vapour_0sal_kPa .* exp(-0.018 * mooring_data.osmotic_coef .* mooring_data.molality);

% Convert to atm
mooring_data.vapour_press_atm = mooring_data.vapour_press_kPa/101.325;

%%
% Henry's law
% As we are unable to know if a particular parcel of water and the gases 
% it contains was in equilibrium with the surface water, or when it was at 
% the surface, we will use standard atmospheric gas fractions to make our 
% calculations. 

% We now calculate the Henry's law constants in um/(kg*atm), where the
% pressure is the partial pressure of the relevant gas. 
% pressure is the partial pressure of the relevant gas. 
% The solubility variables are calculated at measured T,S, assuming a wet standard atmosphere.
% The corresponding partial pressure of the gas (O2 or N2) is calculated from the dry atmosphere mole fractions, 
% and must thus be corrected for  the presence of water vapour, to yield Henry's law constants for the wet atmosphere.


mooring_data.Henry_law_constant_O2 = mooring_data.dox2_sol_umolkg./((1-mooring_data.vapour_press_atm)*constants.mole_fraction_O2);

mooring_data.Henry_law_constant_N2 = mooring_data.N2sol_umolkg./((1-mooring_data.vapour_press_atm)*constants.mole_fraction_N2);

% We now use our own modificaiton of Eqn. 1 from Emerson 2008 to calculate
% the partial pressure of N2 in the water (in atm).

mooring_data.N2_partial_pressure_atm = ((mooring_data.gastension_Pa/(101325)) - mooring_data.vapour_press_atm - (mooring_data.dox2_umolkg./mooring_data.Henry_law_constant_O2))*(constants.mole_fraction_N2/(constants.mole_fraction_N2+constants.mole_fraction_Ar+constants.mole_fraction_CO2));

% We now calculate the timeseries of N2 

mooring_data.N2_umolkg = mooring_data.Henry_law_constant_N2.*mooring_data.N2_partial_pressure_atm;

mooring_data.N2_molm3 = (mooring_data.density_kgm3).*mooring_data.N2_umolkg.*(10^-6);

mooring_data.N2sol_molm3 = (mooring_data.density_kgm3).*mooring_data.N2sol_umolkg.*(10^-6);
