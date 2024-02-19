clear; close all; clc

% look at flux 

path_name = './mat_files/';
var_name  = 'sunfish_data';
load(fullfile(path_name,[var_name,'_oxy_qc.mat']));
dat = sunfish; clear sunfish 
% dat.gridded.pressure_grid=ones(size(dat.gridded.pressure)).*dat.gridded.pressure_grid';

dat.gridded.rho = sw_dens(dat.gridded.salinity,dat.gridded.temperature,dat.gridded.pressure);
dat.gridded.oxygen_adjusted_umolkg = dat.gridded.oxygen_adjusted./(dat.gridded.rho/1000);
dat.gridded.oxygen_raw_umolkg = dat.gridded.oxygen_raw./(dat.gridded.rho/1000);
dat.gridded.sigma_t = sw_dens0(dat.gridded.salinity,dat.gridded.temperature);

%% Gridd data to 1 hr and load era5 data
% RESOLUTION ~2.4 hrs and 1.5 km but we will use 1 hr 
t1 = datenum(2022,01,15); t2 = datenum(2022,05,15); 

id = dat.pressure<10;
FAS_1hr.time  = (t1:1/24:t2)';
FAS_1hr.gl_OXY = mean_interp(dat.dateNum(id),dat.adjusted_oxygen_concentration(id),FAS_1hr.time,0);
FAS_1hr.gl_SST = mean_interp(dat.dateNum(id),dat.temperature(id),FAS_1hr.time,0);
FAS_1hr.gl_SSS = mean_interp(dat.dateNum(id),dat.salinity(id),FAS_1hr.time,0);

load('./mat_files/sunfish_era5_data.mat');
b=load('./mat_files/sunfish_era5_wind_data.mat');
era5.u10 = b.era5.u10;
era5.v10 = b.era5.v10;
era5.sst = b.era5.sst;
clear b

%% FAS 1hr resolution
FAS_1hr.u10    = interp1(era5.time,era5.u10,FAS_1hr.time);
FAS_1hr.v10    = interp1(era5.time,era5.v10,FAS_1hr.time);
FAS_1hr.slp_PA = interp1(era5.time,era5.slp,FAS_1hr.time);    
FAS_1hr.slp_atm= FAS_1hr.slp_PA / 101325;
FAS_1hr.U      = sqrt(FAS_1hr.u10.^2 + FAS_1hr.v10.^2);
FAS_1hr.rh     = interp1(era5.time,era5.rh,FAS_1hr.time);
FAS_1hr.gl_OXY_sol = O2solubility_molal(FAS_1hr.gl_SST,FAS_1hr.gl_SSS,FAS_1hr.slp_PA*0.01);

%% FAS 1d sampled resolution
FAS_1d.time = (t1:1:t2)';
FAS_1d.gl_OXY = interp1(FAS_1hr.time,FAS_1hr.gl_OXY,FAS_1d.time);
FAS_1d.gl_SST = interp1(FAS_1hr.time,FAS_1hr.gl_SST,FAS_1d.time);
FAS_1d.gl_SSS = interp1(FAS_1hr.time,FAS_1hr.gl_SSS,FAS_1d.time);

%% FAS 7d sampled resolution
FAS_7d.time = (t1:7:t2+3.5)';
FAS_7d.gl_OXY = interp1(FAS_1hr.time,FAS_1hr.gl_OXY,FAS_7d.time);
FAS_7d.gl_SST = interp1(FAS_1hr.time,FAS_1hr.gl_SST,FAS_7d.time);
FAS_7d.gl_SSS = interp1(FAS_1hr.time,FAS_1hr.gl_SSS,FAS_7d.time);

%% FAS 14d sampled resolution 
FAS_14d.time = (t1:14:t2+7)';
FAS_14d.gl_OXY = interp1(FAS_1hr.time,FAS_1hr.gl_OXY,FAS_14d.time);
FAS_14d.gl_SST = interp1(FAS_1hr.time,FAS_1hr.gl_SST,FAS_14d.time);
FAS_14d.gl_SSS = interp1(FAS_1hr.time,FAS_1hr.gl_SSS,FAS_14d.time);



%% COMPUTE FLUXES
% Use David Nicholson Gas Toolbox
addpath(genpath('./processing_tools/gas_toolbox-master/'))

FAS_1d.gl_OXY = interp1(FAS_1d.time,FAS_1d.gl_OXY,FAS_1hr.time);
FAS_1d.gl_SST = interp1(FAS_1d.time,FAS_1d.gl_SST,FAS_1hr.time);
FAS_1d.gl_SSS = interp1(FAS_1d.time,FAS_1d.gl_SSS,FAS_1hr.time);
FAS_1d.time = FAS_1hr.time;
FAS_1d.slp_atm = FAS_1hr.slp_atm;
FAS_1d.U = FAS_1hr.U;
FAS_1d.rh = FAS_1hr.rh;

FAS_7d.gl_OXY = interp1(FAS_7d.time,FAS_7d.gl_OXY,FAS_1hr.time);
FAS_7d.gl_SST = interp1(FAS_7d.time,FAS_7d.gl_SST,FAS_1hr.time);
FAS_7d.gl_SSS = interp1(FAS_7d.time,FAS_7d.gl_SSS,FAS_1hr.time);
FAS_7d.time = FAS_1hr.time;
FAS_7d.slp_atm = FAS_1hr.slp_atm;
FAS_7d.U = FAS_1hr.U;
FAS_7d.rh = FAS_1hr.rh;

FAS_14d.gl_OXY = interp1(FAS_14d.time,FAS_14d.gl_OXY,FAS_1hr.time);
FAS_14d.gl_SST = interp1(FAS_14d.time,FAS_14d.gl_SST,FAS_1hr.time);
FAS_14d.gl_SSS = interp1(FAS_14d.time,FAS_14d.gl_SSS,FAS_1hr.time);
FAS_14d.time = FAS_1hr.time;
FAS_14d.slp_atm = FAS_1hr.slp_atm;
FAS_14d.U = FAS_1hr.U;
FAS_14d.rh = FAS_1hr.rh;


% % Stanley et al. 2009
[FAS_1hr.Fd, FAS_1hr.Fc, FAS_1hr.Fp, FAS_1hr.Deq, FAS_1hr.ks] = fas(FAS_1hr.gl_OXY*1e-3,FAS_1hr.U,FAS_1hr.gl_SSS,FAS_1hr.gl_SST,FAS_1hr.slp_atm,'O2','S09',FAS_1hr.rh);
[FAS_1d.Fd, FAS_1d.Fc, FAS_1d.Fp, FAS_1d.Deq, FAS_1d.ks] = fas(FAS_1d.gl_OXY*1e-3,FAS_1d.U,FAS_1d.gl_SSS,FAS_1d.gl_SST,FAS_1d.slp_atm,'O2','S09',FAS_1d.rh);
[FAS_7d.Fd, FAS_7d.Fc, FAS_7d.Fp, FAS_7d.Deq, FAS_7d.ks] = fas(FAS_7d.gl_OXY*1e-3,FAS_7d.U,FAS_7d.gl_SSS,FAS_7d.gl_SST,FAS_7d.slp_atm,'O2','S09',FAS_7d.rh);
[FAS_14d.Fd, FAS_14d.Fc, FAS_14d.Fp, FAS_14d.Deq, FAS_14d.ks] = fas(FAS_14d.gl_OXY*1e-3,FAS_14d.U,FAS_14d.gl_SSS,FAS_14d.gl_SST,FAS_14d.slp_atm,'O2','S09',FAS_14d.rh);

%% INTEGRATE MONTHLY FLUXES

% Define the time period
startDate = datenum('01-Jan-2022');
endDate = datenum('30-May-2022');

% Combine the flux components
FAS_1hr.totalFlux = FAS_1hr.Fd + FAS_1hr.Fc + FAS_1hr.Fp;
FAS_1d.totalFlux = FAS_1d.Fd + FAS_1d.Fc + FAS_1d.Fp;
FAS_7d.totalFlux = FAS_7d.Fd + FAS_7d.Fc + FAS_7d.Fp;
FAS_14d.totalFlux = FAS_14d.Fd + FAS_14d.Fc + FAS_14d.Fp;


% Perform integration for each time series
monthlyFlux_1hr = integrateFlux(FAS_1hr.time,  FAS_1hr.totalFlux, startDate, endDate, t1, t2);
monthlyFlux_1d = integrateFlux(FAS_1d.time,FAS_1d.totalFlux, startDate, endDate, t1, t2);
monthlyFlux_7d = integrateFlux(FAS_7d.time, FAS_7d.totalFlux, startDate, endDate, t1, t2);
monthlyFlux_14d = integrateFlux(FAS_14d.time, FAS_14d.totalFlux, startDate, endDate, t1, t2);



%% Analyze Results
t=tiledlayout(3,1,"TileSpacing","compact",'Padding','compact');

nexttile; hold on
h1=bar(FAS_1hr.time,-(FAS_1hr.Fd+FAS_1hr.Fc+FAS_1hr.Fp)*86400,'stacked','cyan')
h2=plot(FAS_1d.time,-(FAS_1d.Fd+FAS_1d.Fc+FAS_1d.Fp)*86400,'-k')
h3=plot(FAS_7d.time,-(FAS_7d.Fd+FAS_7d.Fc+FAS_7d.Fp)*86400,'-b')
h4=plot(FAS_14d.time,-(FAS_14d.Fd+FAS_14d.Fc+FAS_14d.Fp)*86400,'-r')
hlines(0,'-k')
ylabel('F_{O_2,total} / mol m^2 d^{-1}')
datetick('x','dd-mmm','keeplimits')
ylim([-1.2 0.4])
legend([h1 h2 h3 h4],{'1-hr','1-day','7-day','14-day'},'Location','best');

title('(a)')
formatplot



% Plot the percentage differences


% % Calculate percentage differences relative to the 1-hour resolution
% pct_diff_1d = 100 * (monthlyFlux_1d - monthlyFlux_1hr) ./ monthlyFlux_1hr;
% pct_diff_7d = 100 * (monthlyFlux_7d - monthlyFlux_1hr) ./ monthlyFlux_1hr;
% pct_diff_14d = 100 * (monthlyFlux_14d - monthlyFlux_1hr) ./ monthlyFlux_1hr;
% 
% bar(monthIndices + offsets(1), pct_diff_1d, barWidth, 'FaceColor', 'k');
% bar(monthIndices + offsets(2), pct_diff_7d, barWidth, 'FaceColor', 'b');
% bar(monthIndices + offsets(3), pct_diff_14d, barWidth, 'FaceColor', 'r');
% 
% xlabel('Month');
% ylabel('Percentage Difference (%)');
% title('Percentage Difference Relative to 1hr Resolution');
% legend({'1d vs 1hr', '7d vs 1hr', '14d vs 1hr'}, 'Location', 'best');
% set(gca, 'XTick', monthIndices, 'XTickLabel', months);
% formatplot;
% ylim([-50 50])






t1 = datenum(2022,01,15); t2 = datenum(2022,05,15); 
pg = 0:1:1010; pg = pg(:);
[Xq,Yq]=meshgrid(t1:4/24:t2,pg);


[OXY,idx]=deleteAlmostEmptyColumns(dat.gridded.oxygen_adjusted_umolkg,dat.gridded.pressure_grid);
OXY_intp  = interp2(dat.gridded.timeg(idx),dat.gridded.pressure_grid(1:1010)',...
    OXY(1:1010,:),Xq,Yq);
S  = interp2(dat.gridded.timeg,dat.gridded.pressure_grid(1:1010)',...
    dat.gridded.salinity(1:1010,:),Xq,Yq);
T  = interp2(dat.gridded.timeg,dat.gridded.pressure_grid(1:1010)',...
    dat.gridded.temperature(1:1010,:),Xq,Yq);

tq = dat.gridded.timeg(idx);

timeDiff = diff(tq);
gapStartIndices = find(timeDiff > 1); % 1/4 day threshold
gapEndIndices = gapStartIndices + 1;

% Convert indices to corresponding times
gapStartTimes = tq(gapStartIndices);
gapEndTimes = tq(gapEndIndices);

for i = 1:length(gapStartTimes)
    % Define the time range of the gap
    gapStart = gapStartTimes(i);
    gapEnd = gapEndTimes(i);

    % Find the corresponding indices in the interpolated grid
    gapMask = Xq > gapStart & Xq < gapEnd;

    % Set interpolated data to NaN within the gap
    OXY_intp(gapMask) = NaN;
    T(gapMask)=NaN;
    S(gapMask)=NaN;
end
nprof=length(Xq(1,:));
pmld = zeros(nprof,3);
for ix = 1:nprof
    pmld(ix,:) = mld(Yq(:,ix), T(:,ix), S(:,ix), 'metric', 'threshold', ...
        'tthresh', 0.05, 'dthresh', 0.01,'refpres',25);
end
MLD = pmld(:,3); % density threshold to estimate MLD 
% idt = (Xq(1,:)>datenum(2022,04,15) | Xq(1,:)<datenum(2022,03,01)) & MLD'>999;
% MLD(idt)=NaN;
MLD(MLD>999)=NaN;


id600 = find(abs(Yq(:,1)-600)<5);

nexttile; hold on;
plot(Xq(1,:),-MLD,':k');
plot(Xq(1,:),-MLD,'.r');
ylabel('MLD / m')
title('b')

xlim([FAS_1hr.time(1) FAS_1hr.time(end)])
datetick('x','dd-mmm','keeplimits')
formatplot
yyaxis('right')
plot(Xq(1,:),OXY_intp(id600,:),'-b','LineWidth',1.5);

ylabel('O_2 (595<z<605) / \mumol kg^{-1}')
formatplot
set(gca,'YColor','b')
xlim([FAS_1hr.time(1) FAS_1hr.time(end)])
datetick('x','dd.mmm','keeplimits')



nexttile; hold on

% Define the names of the months covered in your data
months = {'Jan', 'Feb', 'Mar', 'Apr', 'May'};
monthIndices = 1:length(months); % Numeric indices for the months

% Define bar width and offsets
barWidth = 0.15;
offsets = [-1.5, -0.5, 0.5, 1.5] * barWidth;

bar(monthIndices + offsets(1), monthlyFlux_1hr, barWidth, 'FaceColor', 'c');
bar(monthIndices + offsets(2), monthlyFlux_1d, barWidth, 'FaceColor', 'k');
bar(monthIndices + offsets(3), monthlyFlux_7d, barWidth, 'FaceColor', 'b');
bar(monthIndices + offsets(4), monthlyFlux_14d, barWidth, 'FaceColor', 'r');
xlabel('Month');
ylabel('F_{O_2,total} / mol m^2');
title('(c)');
legend({'1hr Resolution', '1d Resolution', '7d Resolution', '14d Resolution'}, 'Location', 'best');
set(gca, 'XTick', monthIndices, 'XTickLabel', months);
formatplot
ylim([0 7])


subtitle(t,'Glider Sunfish Flux Computations')


save_figure(gcf,'./plots/sunfish_O2_flux_comparison',[7.5 7],'.png','300')













function monthlyFlux = integrateFlux(timeSeries, totalFlux, startDate, endDate, t1, t2)
    % Time vector for 1-hour resolution over the entire period
    hourlyTimeVec = (startDate:1/24:endDate)';

    % Interpolate total flux to 1-hour resolution
    id = ~isnan(totalFlux);
    hourlyFlux = interp1(timeSeries(id), totalFlux(id), hourlyTimeVec, 'linear');
    hourlyFlux(isnan(hourlyFlux))=0;

    % Integrate for each month and scale for missing data
    [~,mm,~] = datevec(hourlyTimeVec);
    uniqueMonths = unique(mm);
    monthlyFlux = zeros(length(uniqueMonths), 1);

    TimeVec  = hourlyTimeVec-hourlyTimeVec(1); 
    for i = 1:length(uniqueMonths)
        monthMask = mm == uniqueMonths(i);
        
        % Integrate flux over the month (1-hour resolution, converted to seconds)
        monthlyIntegration = trapz(TimeVec(monthMask), hourlyFlux(monthMask))*86400;

        % Calculate the scaling factor based on available data
        availableDays = sum(hourlyTimeVec(monthMask) >= t1 & hourlyTimeVec(monthMask) <= t2);
        totalDays = sum(monthMask);
        scalingFactor = availableDays / totalDays;

        % Scale the integrated result
        monthlyFlux(i) = monthlyIntegration * scalingFactor;
    end
end














% clear; close all; clc
% 
% % look at fluxes
% 
% path_name = './mat_files/';
% var_name  = 'sunfish_data';
% load(fullfile(path_name,[var_name,'_oxy_qc.mat']));
% dat = sunfish; 
% dat.gridded.pressure_grid=ones(size(dat.gridded.pressure)).*dat.gridded.pressure_grid';
% 
% dat.gridded.rho = sw_dens(dat.gridded.salinity,dat.gridded.temperature,dat.gridded.pressure);
% dat.gridded.oxygen_adjusted_umolkg = dat.gridded.oxygen_adjusted./(dat.gridded.rho/1000);
% dat.gridded.oxygen_raw_umolkg = dat.gridded.oxygen_raw./(dat.gridded.rho/1000);
% dat.gridded.sigma_t = sw_dens0(dat.gridded.salinity,dat.gridded.temperature);
% 
% %% Gridd data to 1 hr and load era5 data
% % RESOLUTION ~2.4 hrs and 1.5 km but we will use 1 hr 
% t1 = datenum(2022,01,1); t2 = datenum(2022,05,30); 
% 
% id = dat.pressure<10;
% FAS_1hr.time  = (t1:1/24:t2)';
% FAS_1hr.gl_OXY = mean_interp(dat.dateNum(id),dat.raw_oxygen_concentration(id),FAS_1hr.time,0);
% FAS_1hr.gl_SST = mean_interp(dat.dateNum(id),dat.temperature(id),FAS_1hr.time,0);
% FAS_1hr.gl_SSS = mean_interp(dat.dateNum(id),dat.salinity(id),FAS_1hr.time,0);
% 
% load('./mat_files/sunfish_era5_data.mat');
% b=load('./mat_files/sunfish_era5_wind_data.mat');
% era5.u10 = b.era5.u10;
% era5.v10 = b.era5.v10;
% era5.sst = b.era5.sst;
% clear b
% 
% %% FAS 1hr resolution
% FAS_1hr.u10    = interp1(era5.time,era5.u10,FAS_1hr.time);
% FAS_1hr.v10    = interp1(era5.time,era5.v10,FAS_1hr.time);
% FAS_1hr.slp_PA = interp1(era5.time,era5.slp,FAS_1hr.time);    
% FAS_1hr.slp_atm= FAS_1hr.slp_PA / 101325;
% FAS_1hr.U      = sqrt(FAS_1hr.u10.^2 + FAS_1hr.v10.^2);
% FAS_1hr.rh     = interp1(era5.time,era5.rh,FAS_1hr.time);
% FAS_1hr.gl_OXY_sol = O2solubility_molal(FAS_1hr.gl_SST,FAS_1hr.gl_SSS,FAS_1hr.slp_PA*0.01);
% 
% 
% %% FAS 1d sampled resolution
% FAS_1d.time = (t1:1:t2)';
% FAS_1d.slp_atm = interp1(FAS_1hr.time,FAS_1hr.slp_atm,FAS_1d.time);
% FAS_1d.U = interp1(FAS_1hr.time,FAS_1hr.U,FAS_1d.time);
% FAS_1d.rh = interp1(FAS_1hr.time,FAS_1hr.rh,FAS_1d.time);
% FAS_1d.gl_OXY = interp1(FAS_1hr.time,FAS_1hr.gl_OXY,FAS_1d.time);
% FAS_1d.gl_SST = interp1(FAS_1hr.time,FAS_1hr.gl_SST,FAS_1d.time);
% FAS_1d.gl_SSS = interp1(FAS_1hr.time,FAS_1hr.gl_SSS,FAS_1d.time);
% 
% %% FAS 7d sampled resolution
% FAS_7d.time = (t1:7:t2+3.5)';
% FAS_7d.slp_atm = interp1(FAS_1hr.time,FAS_1hr.slp_atm,FAS_7d.time);
% FAS_7d.U = interp1(FAS_1hr.time,FAS_1hr.U,FAS_7d.time);
% FAS_7d.rh = interp1(FAS_1hr.time,FAS_1hr.rh,FAS_7d.time);
% FAS_7d.gl_OXY = interp1(FAS_1hr.time,FAS_1hr.gl_OXY,FAS_7d.time);
% FAS_7d.gl_SST = interp1(FAS_1hr.time,FAS_1hr.gl_SST,FAS_7d.time);
% FAS_7d.gl_SSS = interp1(FAS_1hr.time,FAS_1hr.gl_SSS,FAS_7d.time);
% 
% 
% %% FAS 14d sampled resolution 
% FAS_14d.time = (t1:14:t2+7)';
% FAS_14d.slp_atm = interp1(FAS_1hr.time,FAS_1hr.slp_atm,FAS_14d.time);
% FAS_14d.U = interp1(FAS_1hr.time,FAS_1hr.U,FAS_14d.time);
% FAS_14d.rh = interp1(FAS_1hr.time,FAS_1hr.rh,FAS_14d.time);
% FAS_14d.gl_OXY = interp1(FAS_1hr.time,FAS_1hr.gl_OXY,FAS_14d.time);
% FAS_14d.gl_SST = interp1(FAS_1hr.time,FAS_1hr.gl_SST,FAS_14d.time);
% FAS_14d.gl_SSS = interp1(FAS_1hr.time,FAS_1hr.gl_SSS,FAS_14d.time);
% 
% %% COMPUTE FLUXES
% % Use David Nicholson Gas Toolbox
% addpath(genpath('./processing_tools/gas_toolbox-master/'))
% 
% % % Stanley et al. 2009
% [FAS_1hr.Fd, FAS_1hr.Fc, FAS_1hr.Fp, FAS_1hr.Deq, FAS_1hr.ks] = fas(FAS_1hr.gl_OXY*1e-3,FAS_1hr.U,FAS_1hr.gl_SSS,FAS_1hr.gl_SST,FAS_1hr.slp_atm,'O2','S09',FAS_1hr.rh);
% [FAS_1d.Fd, FAS_1d.Fc, FAS_1d.Fp, FAS_1d.Deq, FAS_1d.ks] = fas(FAS_1d.gl_OXY*1e-3,FAS_1d.U,FAS_1d.gl_SSS,FAS_1d.gl_SST,FAS_1d.slp_atm,'O2','S09',FAS_1d.rh);
% [FAS_7d.Fd, FAS_7d.Fc, FAS_7d.Fp, FAS_7d.Deq, FAS_7d.ks] = fas(FAS_7d.gl_OXY*1e-3,FAS_7d.U,FAS_7d.gl_SSS,FAS_7d.gl_SST,FAS_7d.slp_atm,'O2','S09',FAS_7d.rh);
% [FAS_14d.Fd, FAS_14d.Fc, FAS_14d.Fp, FAS_14d.Deq, FAS_14d.ks] = fas(FAS_14d.gl_OXY*1e-3,FAS_14d.U,FAS_14d.gl_SSS,FAS_14d.gl_SST,FAS_14d.slp_atm,'O2','S09',FAS_14d.rh);
% 
% % Fs09 = (Fd+Fc+Fp);
% % Fdiff09 = Fd;
% % Fbubb09 = Fc+Fp;
% 
% 
% %% INTEGRATE MONTHLY FLUXES
% 
% % Define the time period
% startDate = datenum('01-Jan-2022');
% endDate = datenum('30-May-2022');
% 
% % Combine the flux components
% FAS_1hr.totalFlux = FAS_1hr.Fd + FAS_1hr.Fc + FAS_1hr.Fp;
% FAS_1d.totalFlux = FAS_1d.Fd + FAS_1d.Fc + FAS_1d.Fp;
% FAS_7d.totalFlux = FAS_7d.Fd + FAS_7d.Fc + FAS_7d.Fp;
% FAS_14d.totalFlux = FAS_14d.Fd + FAS_14d.Fc + FAS_14d.Fp;
% 
% 
% 
% % Perform integration for each time series
% monthlyFlux_1hr = integrateFlux(FAS_1hr.time,  FAS_1hr.totalFlux, startDate, endDate, t1, t2);
% monthlyFlux_1d = integrateFlux(FAS_1d.time,FAS_1d.totalFlux, startDate, endDate, t1, t2);
% monthlyFlux_7d = integrateFlux(FAS_7d.time, FAS_7d.totalFlux, startDate, endDate, t1, t2);
% monthlyFlux_14d = integrateFlux(FAS_14d.time, FAS_14d.totalFlux, startDate, endDate, t1, t2);
% 
% 
% 
% %% Analyze Results
% t=tiledlayout(3,1,"TileSpacing","compact",'Padding','compact');
% 
% nexttile; hold on
% bar(FAS_1hr.time,-(FAS_1hr.Fd+FAS_1hr.Fc+FAS_1hr.Fp)*86400,'stacked','cyan')
% plot(FAS_1d.time,-(FAS_1d.Fd+FAS_1d.Fc+FAS_1d.Fp)*86400,'-k')
% plot(FAS_7d.time,-(FAS_7d.Fd+FAS_7d.Fc+FAS_7d.Fp)*86400,'*-b')
% plot(FAS_14d.time,-(FAS_14d.Fd+FAS_14d.Fc+FAS_14d.Fp)*86400,'Marker','square','Color','r','LineStyle','-')
% ylabel('F_{O_2,total} / mol m^2 d^{-1}')
% datetick('x','dd-mmm','keeplimits')
% ylim([-1.2 0.4])
% legend('1-hr','1-day','7-day','14-day','Location','best')
% hlines(0,'-k')
% title('(a)')
% formatplot
% 
% nexttile; hold on
% bar(FAS_1hr.time,-(FAS_1hr.Fc+FAS_1hr.Fp)*86400,'stacked','cyan')
% plot(FAS_1d.time,-(FAS_1d.Fc+FAS_1d.Fp)*86400,'-k')
% plot(FAS_7d.time,-(FAS_7d.Fc+FAS_7d.Fp)*86400,'*-b')
% plot(FAS_14d.time,-(FAS_14d.Fc+FAS_14d.Fp)*86400,'Marker','square','Color','r','LineStyle','-')
% ylabel('F_{O_2,bubble} / mol m^2 d^{-1}')
% datetick('x','dd-mmm','keeplimits')
% ylim([-1.2 0.4])
% hlines(0,'-k')
% title('(b)')
% formatplot
% 
% 
% nexttile; hold on
% 
% % Define the names of the months covered in your data
% months = {'Jan', 'Feb', 'Mar', 'Apr', 'May'};
% monthIndices = 1:length(months); % Numeric indices for the months
% 
% % Define bar width and offsets
% barWidth = 0.15;
% offsets = [-1.5, -0.5, 0.5, 1.5] * barWidth;
% 
% bar(monthIndices + offsets(1), monthlyFlux_1hr, barWidth, 'FaceColor', 'c');
% bar(monthIndices + offsets(2), monthlyFlux_1d, barWidth, 'FaceColor', 'k');
% bar(monthIndices + offsets(3), monthlyFlux_7d, barWidth, 'FaceColor', 'b');
% bar(monthIndices + offsets(4), monthlyFlux_14d, barWidth, 'FaceColor', 'r');
% xlabel('Month');
% ylabel('F_{O_2,total} / mol m^2');
% title('(c)');
% legend({'1hr Resolution', '1d Resolution', '7d Resolution', '14d Resolution'}, 'Location', 'best');
% set(gca, 'XTick', monthIndices, 'XTickLabel', months);
% formatplot
% 
% 
% subtitle(t,'Glider Sunfish Flux Computations')
% 
% 
% save_figure(gcf,'./plots/sunfish_O2_flux_comparison',[7.5 7],'.png','300')
% 
% 
% 
% %% Investigate flux dependence on subsampling 
% % Define the range of shifts for the starting point (e.g., 0 to 6 days)
% dshifts = -3:1:3; % 7-day interval
% totalMonthlyFluxShifted = zeros(length(months),length(dshifts));
% 
% % Preallocate array to store total flux for each shift
% % Loop over each shift and calculate total flux
% for i = 1:length(dshifts)
%     % Shift the starting point of the subsampling
%     FAS_7d.time = (t1 + dshifts(i):7:t2)';
%     FAS_7d.slp_atm = interp1(FAS_1hr.time,FAS_1hr.slp_atm,FAS_7d.time);
%     FAS_7d.U = interp1(FAS_1hr.time,FAS_1hr.U,FAS_7d.time);
%     FAS_7d.rh = interp1(FAS_1hr.time,FAS_1hr.rh,FAS_7d.time);
%     FAS_7d.gl_OXY = interp1(FAS_1hr.time,FAS_1hr.gl_OXY,FAS_7d.time);
%     FAS_7d.gl_SST = interp1(FAS_1hr.time,FAS_1hr.gl_SST,FAS_7d.time);
%     FAS_7d.gl_SSS = interp1(FAS_1hr.time,FAS_1hr.gl_SSS,FAS_7d.time);
%     
%     % Calculate total flux using the integrateFlux function
%     [FAS_7d.Fd, FAS_7d.Fc, FAS_7d.Fp, FAS_7d.Deq, FAS_7d.ks] = fas(FAS_7d.gl_OXY*1e-3,FAS_7d.U,FAS_7d.gl_SSS,FAS_7d.gl_SST,FAS_7d.slp_atm,'O2','S09',FAS_7d.rh);
%     FAS_7d.totalFlux = FAS_7d.Fd + FAS_7d.Fc + FAS_7d.Fp;
%     totalMonthlyFluxShifted(:,i) = integrateFlux(FAS_7d.time, FAS_7d.totalFlux, startDate, endDate, t1, t2);
% end
% % Plot total flux as a function of the shift in the starting point
% figure;
% bar(totalMonthlyFluxShifted, 'grouped');
% xlabel('Month');
% ylabel('Total Flux (mol/m²)');
% title('Dependence of Total Monthly Flux on Subsampling Strategy for 7-day Data');
% legend(arrayfun(@(x) sprintf('Shift %d days', x), dshifts, 'UniformOutput', false), 'Location', 'best');
% set(gca, 'XTickLabel', months);
% save_figure(gcf,'./plots/sunfish_flux_subsampling_dependence',[7.5 4],'.png','300')
% 
% 
% % monthlyFlux_7d = nanmean(totalMonthlyFluxShifted,2)
% % dshifts = -7:1:7; % 7-day interval
% % totalMonthlyFluxShifted = zeros(length(months),length(dshifts));
% % % Preallocate array to store total flux for each shift
% % % Loop over each shift and calculate total flux
% % for i = 1:length(dshifts)
% %     FAS_14d.time = (t1 + dshifts(i):7:t2)';
% %     FAS_14d.slp_atm = interp1(FAS_1hr.time,FAS_1hr.slp_atm,FAS_14d.time);
% %     FAS_14d.U = interp1(FAS_1hr.time,FAS_1hr.U,FAS_14d.time);
% %     FAS_14d.rh = interp1(FAS_1hr.time,FAS_1hr.rh,FAS_14d.time);
% %     FAS_14d.gl_OXY = interp1(FAS_1hr.time,FAS_1hr.gl_OXY,FAS_14d.time);
% %     FAS_14d.gl_SST = interp1(FAS_1hr.time,FAS_1hr.gl_SST,FAS_14d.time);
% %     FAS_14d.gl_SSS = interp1(FAS_1hr.time,FAS_1hr.gl_SSS,FAS_14d.time);
% %     
% %     % Calculate total flux using the integrateFlux function
% %     [FAS_14d.Fd, FAS_14d.Fc, FAS_14d.Fp, FAS_14d.Deq, FAS_14d.ks] = fas(FAS_14d.gl_OXY*1e-3,FAS_14d.U,FAS_14d.gl_SSS,FAS_14d.gl_SST,FAS_14d.slp_atm,'O2','S09',FAS_14d.rh);
% %     FAS_14d.totalFlux = FAS_14d.Fd + FAS_14d.Fc + FAS_14d.Fp;
% %     totalMonthlyFluxShifted(:,i) = integrateFlux(FAS_14d.time, FAS_14d.totalFlux, startDate, endDate, t1, t2);
% % end
% % monthlyFlux_14d = nanmean(totalMonthlyFluxShifted,2)
% % 
% % 
% % 
% % figure(); hold on
% % % Define the names of the months covered in your data
% % months = {'Jan', 'Feb', 'Mar', 'Apr', 'May'};
% % monthIndices = 1:length(months); % Numeric indices for the months
% % 
% % % Define bar width and offsets
% % barWidth = 0.15;
% % offsets = [-1.5, -0.5, 0.5, 1.5] * barWidth;
% % 
% % bar(monthIndices + offsets(1), monthlyFlux_1hr, barWidth, 'FaceColor', 'c');
% % bar(monthIndices + offsets(2), monthlyFlux_1d, barWidth, 'FaceColor', 'k');
% % bar(monthIndices + offsets(3), monthlyFlux_7d, barWidth, 'FaceColor', 'b');
% % bar(monthIndices + offsets(4), monthlyFlux_14d, barWidth, 'FaceColor', 'r');
% % xlabel('Month');
% % ylabel('F_{O_2,total} / mol m^2');
% % title('(c)');
% % legend({'1hr Resolution', '1d Resolution', '7d Resolution', '14d Resolution'}, 'Location', 'best');
% % set(gca, 'XTick', monthIndices, 'XTickLabel', months);
% % formatplot
% % save_figure(gcf,'./plots/sunfish_O2_monthly_flux',[7.5 7],'.png','300')
% % 
% 
% 
% 
% 
% 
% function monthlyFlux = integrateFlux(timeSeries, totalFlux, startDate, endDate, t1, t2)
%     % Time vector for 1-hour resolution over the entire period
%     hourlyTimeVec = (startDate:1/24:endDate)';
% 
%     % Interpolate total flux to 1-hour resolution
%     id = ~isnan(totalFlux);
%     hourlyFlux = interp1(timeSeries(id), totalFlux(id), hourlyTimeVec, 'linear');
%     hourlyFlux(isnan(hourlyFlux))=0;
% 
%     % Integrate for each month and scale for missing data
%     [~,mm,~] = datevec(hourlyTimeVec);
%     uniqueMonths = unique(mm);
%     monthlyFlux = zeros(length(uniqueMonths), 1);
% 
%     TimeVec  = hourlyTimeVec-hourlyTimeVec(1); 
%     for i = 1:length(uniqueMonths)
%         monthMask = mm == uniqueMonths(i);
%         
%         % Integrate flux over the month (1-hour resolution, converted to seconds)
%         monthlyIntegration = trapz(TimeVec(monthMask), hourlyFlux(monthMask))*86400;
% 
%         % Calculate the scaling factor based on available data
%         availableDays = sum(hourlyTimeVec(monthMask) >= t1 & hourlyTimeVec(monthMask) <= t2);
%         totalDays = sum(monthMask);
%         scalingFactor = availableDays / totalDays;
% 
%         % Scale the integrated result
%         monthlyFlux(i) = monthlyIntegration * scalingFactor;
%     end
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
