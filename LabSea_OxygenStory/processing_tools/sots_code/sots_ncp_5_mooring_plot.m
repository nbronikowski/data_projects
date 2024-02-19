% Plots the data contained in mooring_data

%% Assign atmospheric pressure (allows solubility plotting)

% If a timeseries of atmospheric pressure is not available in mooring_data, 
% we use the user defined steady choice

if ~isfield(mooring_data,'atmosphericpress_Pa') || atmospheric_pressure_manual_override

    mooring_data.atmosphericpress_Pa = (constants.atm_in_Pa*atmospheric_pressure_choice) * ones(size(mooring_data.time));

    disp(['Atmospheric pressure: Constant user choice of ',num2str(atmospheric_pressure_choice),'atm used.'])

else

    disp(strcat('Atmospheric pressure: Timeseries available from mooring used'))

end

%% Grouped plots

% Produces a plot of temperature and salinity
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
yyaxis left
plot(mooring_data.time,mooring_data.temp_C)
ylabel('Temp °C')
yyaxis right
plot(mooring_data.time,mooring_data.psal_PSU)
ylabel(' Salinity PSU')
xlim([mooring_data.time(1) mooring_data.time(end)])
datetick(gca,'KeepLimits')
title('Temperature and salinity')
grid on

% Produces a plot of measured oxygen and oxygen solubility
subplot(2,2,2)
yyaxis left
plot(mooring_data.time,mooring_data.dox2_umolkg)
ylabel('Measured O2 umol/kg')
yyaxis right
plot(mooring_data.time,(mooring_data.atmosphericpress_Pa/constants.atm_in_Pa).*mooring_data.dox2_sol_umolkg)
title('Oxygen')
ylabel('O2 solubility umol/kg')
xlim([mooring_data.time(1) mooring_data.time(end)])
datetick(gca,'KeepLimits')
grid on

% Produces a plot of gas tension and water density
subplot(2,2,3)
yyaxis left
plot(mooring_data.time,mooring_data.gastension_Pa)
ylabel('Gas tension Pa')
yyaxis right
plot(mooring_data.time,mooring_data.density_kgm3)
title('Gas tension and density')
ylabel('Density kg/m3')
xlim([mooring_data.time(1) mooring_data.time(end)])
datetick(gca,'KeepLimits')
grid on

% Produces a plot of mixed layer depth and windspeed
subplot(2,2,4)
yyaxis left
plot(mooring_data.time,mooring_data.mld_m)
ylabel('Mixed layer depth m')
set(gca,'YDir','reverse')
yyaxis right
plot(mooring_data.time,mooring_data.windspeed_ms)
title('Mixed layer depth and windspeed')
ylabel('10m windspeed m/s')
xlim([mooring_data.time(1) mooring_data.time(end)])
datetick(gca,'KeepLimits')
grid on


