clear; close all; clc

path_name = './mat_files/';
var_name  = 'pearldiver_data';
load(fullfile(path_name,[var_name,'_oxy_qc.mat']));
dat = pearldiver; clear pearldiver;

%% Gridd data to 1 hr and load era5 data
% RESOLUTION ~2.4 hrs and 1.5 km but we will use 1 hr to keep thing simple
t1 = datenum(2020,01,01); t2 = datenum(2020,05,30); 
pg = 0:1:1030; pg = pg(:);
[Xq,Yq]=meshgrid(t1:1/24:t2,pg);

OXY  = interp2(dat.gridded.timeg,dat.gridded.pressure_grid,dat.gridded.oxygen_adjusted,Xq,Yq);
SALT = interp2(dat.gridded.timeg,dat.gridded.pressure_grid,dat.gridded.salinity,Xq,Yq);
TEMP = interp2(dat.gridded.timeg,dat.gridded.pressure_grid,dat.gridded.temperature,Xq,Yq);

%% Reintroduce gaps
timeDiff = diff(dat.gridded.timeg);
gapStartIndices = find(timeDiff > 0.25); % 1/4 day threshold
gapEndIndices = gapStartIndices + 1;

% Convert indices to corresponding times
gapStartTimes = dat.gridded.timeg(gapStartIndices);
gapEndTimes = dat.gridded.timeg(gapEndIndices);

for i = 1:length(gapStartTimes)
    % Define the time range of the gap
    gapStart = gapStartTimes(i);
    gapEnd = gapEndTimes(i);

    % Find the corresponding indices in the interpolated grid
    gapMask = Xq > gapStart & Xq < gapEnd;

    % Set interpolated data to NaN within the gap
    SALT(gapMask) = NaN;
    TEMP(gapMask) = NaN;
end

%% Compute MLD 
nprof=length(Xq(1,:));
pmld = zeros(nprof,3);
for ix = 1:nprof
    pmld(ix,:) = mld(Yq(:,ix), TEMP(:,ix), SALT(:,ix), 'metric', 'threshold', ...
        'tthresh', 0.05, 'dthresh', 0.01);
end
MLD = pmld(:,3); % density threshold to estimate MLD 
idt = (Xq(1,:)>datenum(2020,04,15) | Xq(1,:)<datenum(2020,03,01)) & MLD'>1000;
MLD(idt)=NaN;

% figure()
% imagescn(Xq,-Yq,OXY);
% hold on 
% plot(Xq(1,:),-MLD,'k');cb=colorbar;
% ylabel(cb,'oxygen concentration / \mumol L^{-1}')
% ylabel('depth (m)')
% datetick('x')


%% Compute dO2 / dt 
dO2 = diff(OXY*0.001, 1, 2); % micro mol / L to mol / m3
dt = nanmedian(diff(Xq(1,:)))*86400; % s / hr 
dO2dt = dO2 / dt; % we want rate of change of Oxygen in mol m^-3 s^-1
dO2dt = [dO2dt, NaN(size(dO2dt, 1), 1)];

% dO2dt=medfilt2(dO2dt,[3 10],'symmetric');

figure()

subplot(211)
imagescn(Xq,-Yq,dO2dt);
hold on 
plot(Xq(1,:),-MLD,'k');cb=colorbar;
ylabel(cb,'oxygen concentration / \mumol L^{-1}')
ylabel('depth (m)')
datetick('x')

subplot(212);
imagescn(Xq,-Yq,medfilt2(dO2dt,[20 48],'symmetric'));
hold on 
plot(Xq(1,:),-MLD,'k');cb=colorbar;
ylabel(cb,'oxygen concentration / \mumol L^{-1}')
ylabel('depth (m)')
datetick('x')


dO2dt_filt = medfilt2(dO2dt,[20 48],'symmetric');
dO2dt_filt(gapMask)=NaN;
OXY(gapMask)=NaN;


%% Compute average rates below and above MLD 
O2rateMLD = NaN*MLD;
O2rateBMLD = NaN*MLD;
for i = 1:nprof
    nnan = find(~isnan(dO2dt_filt(:,i)));
    if length(nnan)>3 && ~isnan(MLD(i))
        idMLD = find(Yq(:,1)<(MLD(i)-5));
        idBMLD = find(Yq(:,1)>(MLD(i)+5));
        O2rateMLD(i)=nanmean(dO2dt_filt(idMLD,i)); % rate in mixed layer
        O2rateBMLD(i)=nanmean(dO2dt_filt(idBMLD,i)); % rate below mixed layer
    end
end

figure(); hold on
plot(Xq(1,:),O2rateMLD*86400,'b');
plot(Xq(1,:),O2rateBMLD*86400,'r');
datetick('x');
ylabel('$\frac{dO_2}{dt}$ (mol m$^{-3}$ d$^{-1}$)','Interpreter','latex');


%% OLD METHOD
% id =  dat.dateNum>t1 & dat.dateNum<t2;
% subs.time = (t1:1/1440:t2)';
% 
% subs.pres = mean_interp(dat.dateNum(id),dat.pressure(id),subs.time);
% subs.temp = mean_interp(dat.dateNum(id),dat.temperature(id),subs.time);
% subs.salt = mean_interp(dat.dateNum(id),dat.salinity(id),subs.time);
% subs.oxy  = mean_interp(dat.dateNum(id),dat.adjusted_oxygen_concentration(id),subs.time);
% 
% x2 = scale_var(subs.time,1/64);
% [temp_pg,ux]  = pgrid_columns(x2,subs.pres,subs.temp,pg);
% [salt_pg,~]  = pgrid_columns(x2,subs.pres,subs.salt,pg);
% [oxy_pg,~]  = pgrid_columns(x2,subs.pres,subs.oxy,pg);
% 
% % interpolate to x grid + smooth 
% Nr = 1; Nc = 1;
% OXY =intp_x_dim(ux,oxy_pg,Xq,Yq,Nr,Nc);
% SALT=intp_x_dim(ux,salt_pg,Xq,Yq,Nr,Nc);
% TEMP=intp_x_dim(ux,temp_pg,Xq,Yq,Nr,Nc);
% 
% subplot(211);
% imagescn(OXY)
% 
% subplot(212);
% imagescn(medfilt2(OXY,[2 6],'symmetric'))
% 
% OXY = medfilt2(OXY,[2 6],'symmetric');
