% We calculate schmidt numbers for water of salinity=35 and temperatures 
% range 0-30C (Wanninkhof 92)
mooring_data.Sc_N2_W92 = constants.Sc_A_N2_W92 - constants.Sc_B_N2_W92*mooring_data.temp_C + constants.Sc_C_N2_W92*(mooring_data.temp_C.^2) - constants.Sc_D_N2_W92*(mooring_data.temp_C.^3);

mooring_data.Sc_O2_W92 = constants.Sc_A_O2_W92 - constants.Sc_B_O2_W92*mooring_data.temp_C + constants.Sc_C_O2_W92*(mooring_data.temp_C.^2) - constants.Sc_D_O2_W92*(mooring_data.temp_C.^3);

% We assign the Schmidt values
mooring_data.Sc_N2 = mooring_data.Sc_N2_W92;

mooring_data.Sc_O2 = mooring_data.Sc_O2_W92;

% We calculate the gas transfer velocity using Ho (2006)
mooring_data.Gc_O2_cmhr = (0.266.*mooring_data.windspeed_ms.^2).*((mooring_data.Sc_O2./constants.schmidt_CO2_20C_35PSU).^-0.5);

mooring_data.Gc_N2_cmhr = (0.266.*mooring_data.windspeed_ms.^2).*((mooring_data.Sc_N2./constants.schmidt_CO2_20C_35PSU).^-0.5);
