clear; close all; clc

%% COMPUTE FLUX + INVENT
path_name = './mat_files/';
var_name  = 'pearldiver_data';
load(fullfile(path_name,[var_name,'_oxy_qc.mat']));
dat = pearldiver; 
dat.gridded.pressure_grid=ones(size(dat.gridded.pressure)).*dat.gridded.pressure_grid';

dat.gridded.rho = sw_dens(dat.gridded.salinity,dat.gridded.temperature,dat.gridded.pressure);
dat.gridded.oxygen_adjusted_umolkg = dat.gridded.oxygen_adjusted./(dat.gridded.rho/1000);
dat.gridded.oxygen_raw_umolkg = dat.gridded.oxygen_raw./(dat.gridded.rho/1000);
dat.gridded.sigma_t = sw_dens0(dat.gridded.salinity,dat.gridded.temperature);

%% Gridd data to 1 hr and load era5 data
t1 = datenum(2020,01,1); t2 = datenum(2020,05,30); 
id = dat.pressure<10;
FAS_1hr.time  = (t1:1/24:t2)';
FAS_1hr.gl_OXY = mean_interp(dat.dateNum(id),dat.adjusted_oxygen_concentration(id),FAS_1hr.time,0);
FAS_1hr.gl_SST = mean_interp(dat.dateNum(id),dat.temperature(id),FAS_1hr.time,0);
FAS_1hr.gl_SSS = mean_interp(dat.dateNum(id),dat.salinity(id),FAS_1hr.time,0);

load('./mat_files/pearldiver_era5_data.mat');
b=load('./mat_files/pearldiver_era5_wind_data.mat');
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


%% COMPUTE FLUXES
% Use David Nicholson Gas Toolbox
addpath(genpath('./processing_tools/gas_toolbox-master/'))

% % Stanley et al. 2009
[FAS_1hr.Fd, FAS_1hr.Fc, FAS_1hr.Fp, FAS_1hr.Deq, FAS_1hr.ks] = ...
    fas(FAS_1hr.gl_OXY*1e-3,FAS_1hr.U,FAS_1hr.gl_SSS,FAS_1hr.gl_SST,...
    FAS_1hr.slp_atm,'O2','S09',FAS_1hr.rh);

%% INTEGRATE MONTHLY FLUXES

% Define the time period
startDate = datenum('01-January-2020');
endDate = datenum('30-May-2020');

% Combine the flux components
FAS_1hr.totalFlux = FAS_1hr.Fd + FAS_1hr.Fc + FAS_1hr.Fp;

% Perform integration for each time series
monthlyFlux_1hr = integrateFlux(FAS_1hr.time,  FAS_1hr.totalFlux, startDate, endDate, t1, t2);

%% COLUMN INVENT
t1 = datenum(2020,01,1); t2 = datenum(2020,05,30); 
pg = 0:1:1000; pg = pg(:);
[Xq,Yq]=meshgrid(t1:2/24:t2,pg);

[OXY,idx]=deleteAlmostEmptyColumns(dat.gridded.oxygen_adjusted,dat.gridded.pressure_grid);
O2  = interp2(dat.gridded.timeg(idx),dat.gridded.pressure_grid(1:1010)',...
    OXY(1:1010,:),Xq,Yq);
S  = interp2(dat.gridded.timeg,dat.gridded.pressure_grid(1:1010)',...
    dat.gridded.salinity(1:1010,:),Xq,Yq);
T  = interp2(dat.gridded.timeg,dat.gridded.pressure_grid(1:1010)',...
    dat.gridded.temperature(1:1010,:),Xq,Yq);


O2_sm = movmean(O2,5*12,2,"omitnan","Endpoints","fill");

tq = dat.gridded.timeg(idx);
timeDiff = diff(tq);
gapStartIndices = find(timeDiff > 1); % 1/4 day threshold
gapEndIndices = gapStartIndices + 1;
gapStartTimes = tq(gapStartIndices);
gapEndTimes = tq(gapEndIndices);
for i = 1:length(gapStartTimes)
    gapStart = gapStartTimes(i);
    gapEnd = gapEndTimes(i);
    gapMask = Xq > gapStart & Xq < gapEnd;
    O2_sm(gapMask) = NaN;
    O2(gapMask) = NaN;
    O2_sm_dt(gapMask) = NaN;
    T(gapMask)=NaN;
    S(gapMask)=NaN;
end
nprof=length(Xq(1,:));
pmld = zeros(nprof,3);
for ix = 1:nprof
    pmld(ix,:) = mld(Yq(:,ix), T(:,ix), S(:,ix), 'metric', 'threshold', ...
        'tthresh', 0.05, 'dthresh', 0.01);
end
MLD = pmld(:,3); clear pmld
idt = (Xq(1,:)>datenum(2020,04,15) | Xq(1,:)<datenum(2020,03,01)) & MLD'>999;
MLD(idt)=NaN;

% Integrate O2 gradient
O2_sm_dz_bmld= NaN*Xq(1,:);
O2_sm_dz_mld= NaN*Xq(1,:);
for i = 1:length(Xq(1,:))
    v = O2(:,i);
    z = Yq(:,i);

    idnan = isnan(v);
    zn = z(~idnan);
    z_start = min(zn); z_end = max(zn);
    z_mld  = MLD(i);
    if isnan(z_mld) && Xq(1,i)<datenum(2020,03,01)
        if i>2 && ~isnan(MLD(i-1))
            z_mld = MLD(i-1);
        elseif ~isnan(MLD(i+1))
            z_mld = MLD(i+1);
        end
    end
    
   if isnan(z_mld) && Xq(1,i)>datenum(2020,03,01)
        if ~isnan(MLD(i-1)) && i<length(Xq(1,:))
            z_mld = MLD(i-1);
        elseif ~isnan(MLD(i+1)) && i<length(Xq(1,:))-1
            z_mld = MLD(i+1);
        end
    end

    if (length(v(~idnan))>2) && (z_end-z_start>200) && ...
            z_start<20 && z_end>900 && ~isnan(z_mld)
        if z_end<z_mld; z_end=z_mld+10; end
        if z_start>z_mld; z_start=z_mld-10; end
        zi_mld   = z_start:1:z_mld;
        vi_mld = interp1(zn,v(~idnan),zi_mld)*1e-3;
        zi_bmld = z_mld:1:z_end;
        vi_bmld = interp1(zn,v(~idnan),zi_bmld)*1e-3;
        
%         vi_mld = [0,diff(vi_mld)];
%         vi_bmld = [0,diff(vi_bmld)];

        if length(vi_bmld(~isnan(vi_bmld)))>2
            O2_sm_dz_bmld(i) = trapz(zi_bmld,vi_bmld);
        else
            O2_sm_dz_bmld(i) = NaN;
        end
        if length(vi_mld(~isnan(vi_mld)))>2
            O2_sm_dz_mld(i) = trapz(zi_mld,vi_mld);
        else
            O2_sm_dz_mld(i) = NaN;
        end
    else
        O2_sm_dz_bmld(i) = NaN;
        O2_sm_dz_mld(i) = NaN;
    end
end
id = ~isnan(O2_sm_dz_bmld);
FAS_1hr.time_3d = t1:1/24:t2;
O2_sm_dz_bmld = mean_interp(Xq(1,:),O2_sm_dz_bmld,FAS_1hr.time_3d);
O2_sm_dz_mld = mean_interp(Xq(1,:),O2_sm_dz_mld,FAS_1hr.time_3d);
mld_i = mean_interp(Xq(1,:),MLD,FAS_1hr.time_3d);


idt = MLD>999;
gap_times = Xq(1,idt);

idx = 1:length(FAS_1hr.time_3d);
idx_gap = interp1(FAS_1hr.time_3d,idx,gap_times,'nearest');
O2_sm_dz_bmld(idx_gap)=NaN;
O2_sm_dz_mld(idx_gap)=NaN;
mld_i(idx_gap)=NaN;

O2_cont_MLD = O2_sm_dz_mld;
O2_cont_bMLD = O2_sm_dz_bmld;
O2_cont_MLD_3d = medfilt1(O2_cont_MLD,3*24);
O2_cont_bMLD_3d = medfilt1(O2_cont_bMLD,3*24,'includenan','truncate');

%% PLOT Results
t=tiledlayout(2,1,"TileSpacing","compact",'Padding','compact');

nexttile

yyaxis('right'); hold on
months = datenum(2020,1:5,15);
bar(months, monthlyFlux_1hr,0.15,'FaceColor',rgb('orange'),'FaceAlpha',0.5);
xlim([t1 t2])
ylabel('O_{2,Total} / mol m^2');
text(months(3)+15,7,['Feb - Apr Total = ',...
    num2str(sum(monthlyFlux_1hr(2:4))),' mol m^{2}']);
ylim([-10 10])
title('(a) Daily Wind-driven and Integrated Monthly Flux');
formatplot
set(gca,'YColor',rgb('orange'));


yyaxis('left'); hold on
b1=bar(FAS_1hr.time,-FAS_1hr.totalFlux*86400,'stacked','b');
b2=plot(FAS_1hr.time,-FAS_1hr.totalFlux*86400,'-k');
hlines(0,'-k')
ylabel('F_{O_{2}} / mol m^2 d^{-1}')
xlim([t1 t2])
datetick('x','dd-mmm','keeplimits')

ylim([-1 1])
formatplot
set(gca,'YColor',rgb('blue'));

nexttile; hold on;
yyaxis('right'); hold on
plot(FAS_1hr.time_3d,mld_i,'r');
plot(FAS_1hr.time_3d,medfilt1(mld_i,3*24,'includenan','truncate'),'-r','LineWidth',2);
ylabel('MLD / m')
xlim([FAS_1hr.time(1) FAS_1hr.time(end)])
datetick('x','dd-mmm','keeplimits')
formatplot
set(gca,'YColor','r')


yyaxis('left'); hold on
h1=plot(FAS_1hr.time_3d,O2_cont_bMLD,'-b');
h2=plot(FAS_1hr.time_3d,O2_cont_bMLD_3d,'-b','LineWidth',1.5);

idt1 = FAS_1hr.time_3d>datenum(2020,01,1) & FAS_1hr.time_3d<datenum(2020,01,20);
idt2 = FAS_1hr.time_3d>datenum(2020,04,01) & FAS_1hr.time_3d<datenum(2020,05,30);

hlines(nanmean(O2_cont_bMLD_3d(idt1)))
hlines(nanmean(O2_cont_bMLD_3d(idt2)))
delta_O2cont = nanmean(O2_cont_bMLD_3d(idt2))-nanmean(O2_cont_bMLD_3d(idt1));
text(months(4),200,['\Delta O_2 = ',...
    num2str(delta_O2cont),' mol m^{2}']);
ylabel('O_2 Content / mol m^{2}')
formatplot
set(gca,'YColor','b')
xlim([t1 t2])
datetick('x','dd.mmm','keeplimits')
title('(b) Oxygen Column Inventory Below Mixed Layer')

legend([h1 h2],{'1-hr','3-d'},'Location','East')

subtitle(t,'Glider Pearldiver Flux Computations')


save_figure(gcf,'./plots/pearldiver_O2_flux_comparison',[7.5 5],'.png','300')












function monthlyFlux = integrateFlux(x, y, startDate, endDate, t1, t2)
    % Time vector for 1-hour resolution over the entire period
    xi = (startDate:1/24:endDate)';

    % Interpolate total flux to 1-hour resolution
    id = ~isnan(y); 
    yi = interp1(x(id), y(id), xi, 'linear'); 
    yi(isnan(yi))=0;

    % Integrate for each month and scale for missing data
    [~,mm,~] = datevec(xi);
    uniqueMonths = unique(mm);
    monthlyFlux = zeros(length(uniqueMonths), 1);

    TimeVec  = xi-xi(1); 
    for i = 1:length(uniqueMonths)
        monthMask = mm == uniqueMonths(i);
        
        % Integrate flux over the month (1-hour resolution, converted to seconds)
        monthlyIntegration = trapz(TimeVec(monthMask), yi(monthMask))*86400;

        % Calculate the scaling factor based on available data
        availableDays = sum(xi(monthMask) >= t1 & xi(monthMask) <= t2);
        totalDays = sum(monthMask);
        scalingFactor = availableDays / totalDays;

        % Scale the integrated result
        monthlyFlux(i) = monthlyIntegration * scalingFactor;
    end
end












