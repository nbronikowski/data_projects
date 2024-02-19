% Plots oxygen saturation with temperature and gas tension

%%
% Makes a fullscreen figure
figure('units','normalized','outerposition',[0 0 1 1])

% Plots temperature on the right axis
subplot(2,1,1)
yyaxis right
plot(mooring_data.time,mooring_data.temp_C);
ylabel('Temp ^oC');

% Plots oxygen saturation on the left axis
yyaxis left
plot(mooring_data.time,mooring_data.dox2_sat,'LineWidth',2);
ylabel('O2 saturation');
xlim([mooring_data.time(1) mooring_data.time(end)]);
datetick('x','mmm','KeepLimits');
title('Temperature and oxygen saturation')
grid on;

% Plots gas tension on right axis
subplot(2,1,2)
yyaxis right
plot(mooring_data.time,mooring_data.gastension_Pa);
ylabel('Gas tension Pa');

% Plots oxygen saturation on the left axis
yyaxis left
plot(mooring_data.time,mooring_data.dox2_sat,'LineWidth',2);
ylabel('O2 saturation');
xlim([mooring_data.time(1) mooring_data.time(end)]);
datetick('x','mmm','KeepLimits');
title('Gas tension and oxygen saturation')
grid on;