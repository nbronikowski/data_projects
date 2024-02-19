% We calculate schmidt numbers for water of salinity=35 and temperatures 
% range 0-30C (Wanninkhof 92)
mooring_data.Sc_N2_W92 = constants.Sc_A_N2_W92 - constants.Sc_B_N2_W92*mooring_data.temp_C + constants.Sc_C_N2_W92*(mooring_data.temp_C.^2) - constants.Sc_D_N2_W92*(mooring_data.temp_C.^3);

mooring_data.Sc_O2_W92 = constants.Sc_A_O2_W92 - constants.Sc_B_O2_W92*mooring_data.temp_C + constants.Sc_C_O2_W92*(mooring_data.temp_C.^2) - constants.Sc_D_O2_W92*(mooring_data.temp_C.^3);

% We assign the Schmidt values
mooring_data.Sc_N2 = mooring_data.Sc_N2_W92;

mooring_data.Sc_O2 = mooring_data.Sc_O2_W92;


% We calculate the gas transfer velocity using Wanninkhof (2002)
for i = 1:length(mooring_data.windspeed_ms)

    if mooring_data.windspeed_ms(i)<=6

        mooring_data.Gc_O2_cmhr(i) = (0.31.*mooring_data.windspeed_ms(i).^2).*((mooring_data.Sc_O2(i)./constants.schmidt_CO2_20C_35PSU).^-0.5);

        mooring_data.Gc_N2_cmhr(i) = (0.31.*mooring_data.windspeed_ms(i).^2).*((mooring_data.Sc_N2(i)./constants.schmidt_CO2_20C_35PSU).^-0.5);

    else

        mooring_data.Gc_O2_cmhr(i) = (0.39.*mooring_data.windspeed_ms(i).^2).*((mooring_data.Sc_O2(i)./constants.schmidt_CO2_20C_35PSU).^-0.5);

        mooring_data.Gc_N2_cmhr(i) = (0.39.*mooring_data.windspeed_ms(i).^2).*((mooring_data.Sc_N2(i)./constants.schmidt_CO2_20C_35PSU).^-0.5);

    end

end
