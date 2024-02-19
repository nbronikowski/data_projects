clear; close all; clc

% read in the data and extract time series
% plot those timeseries below each other
% add T-S for select profiles and show Oxygen Properties

path_name = './mat_files/';
var_name  = 'sunfish_data';
load(fullfile(path_name,[var_name,'_oxy_qc.mat']));
dat = sunfish; 
dat.gridded.pressure_grid=ones(size(dat.gridded.pressure)).*dat.gridded.pressure_grid';

dat.gridded.rho = sw_dens(dat.gridded.salinity,dat.gridded.temperature,dat.gridded.pressure);
dat.gridded.oxygen_adjusted_umolkg = dat.gridded.oxygen_adjusted./(dat.gridded.rho/1000);
dat.gridded.oxygen_raw_umolkg = dat.gridded.oxygen_raw./(dat.gridded.rho/1000);
dat.gridded.sigma_t = sw_dens0(dat.gridded.salinity,dat.gridded.temperature);

%% Extract points at 2 quadrants
x1=-51.024; y1=53.7008; % away from boundary
x2=-50.261; y2=54.0585; % closer to boundary
t1 = datenum(2022,01,15);
t2 = datenum(2022,05,15);

dkm1=lldistkm([dat.gridded.lat(:),dat.gridded.lon(:)],...
            [y1,x1]);
dkm2=lldistkm([dat.gridded.lat(:),dat.gridded.lon(:)],...
            [y2,x2]);

idx1 = find(abs(dkm1)<25 & dat.gridded.timeg(:)>t1 & dat.gridded.timeg(:)<t2);
idx2 = find(abs(dkm2)<25 & dat.gridded.timeg(:)>t1 & dat.gridded.timeg(:)<t2);

%% Create 2 structs of timeseries
flist = fieldnames(dat.gridded);
for i = 1:length(flist)
    dat1.(flist{i}) = dat.gridded.(flist{i})(:,idx1);
    dat2.(flist{i}) = dat.gridded.(flist{i})(:,idx2);
end


%% Plot these on a regular time grid





figure()
subplot(211); 
plot(dat1.timeg,dat1.oxygen_adjusted_umolkg(10,:),'.r'); hold on
plot(dat2.timeg,dat2.oxygen_adjusted_umolkg(10,:),'.b');
xlim([t1 t2])
datetick('x','keeplimits')


subplot(212);
plot(dat1.timeg,dat1.oxygen_adjusted_umolkg(800,:),'.r'); hold on
plot(dat2.timeg,dat2.oxygen_adjusted_umolkg(800,:),'.b');
xlim([t1 t2])
datetick('x','keeplimits')


figure()
subplot(121);
scatter(dat1.salinity(:),dat1.temperature(:),20,dat1.oxygen_adjusted_umolkg(:),'filled')
caxis([280 320])
colorbar

subplot(122);
scatter(dat2.salinity(:),dat2.temperature(:),20,dat2.oxygen_adjusted_umolkg(:),'filled')
caxis([280 320])
colorbar
