clear; clc; close all

flist = dir('*.nc');
flist(1)=[];
ttime = double(ncread(flist.name,'time'))/24+datenum(1900,1,1);
lon  = double(ncread(flist.name,'longitude'));
lat  = double(ncread(flist.name,'latitude'));
slp  = double(ncread(flist.name,'msl')); % mean sea level pressure Pa
sp   = double(ncread(flist.name,'sp'));  % surface pressure Pa
t2m  = double(ncread(flist.name,'t2m')); % 2 m air temp K
d2m  = double(ncread(flist.name,'d2m')); % 2 m dew point temp K

%% compute relative humidity
% a1=611.21;a3=17.502;a4=32.19;T0=273.16;
% Esat=a1*exp(a3*(d2m-T0)./(d2m-a4));% 
% Rdry=287.0597;  Rvap=461.5250; 
% qsat=(Rdry/Rvap)*E/(ps-((1-Rdry/Rvap)*E))

d2mC = d2m-273.15;
t2mC = t2m-273.15;
e = 0.6108 * exp((17.27 * d2mC) ./ (d2mC + 237.3)); 
es= 0.6108 * exp((17.27 * t2mC) ./ (t2mC + 237.3));

RH = e./es;

% extract sunfish data
path_name = './../mat_files/'
var_name  = 'sunfish_data'
load(fullfile(path_name,[var_name,'_clean.mat']))
dat = sunfish;  clear sunfish
dat.time = datenum(dat.time);

t1 = datenum(dat.time(1));
t2 = datenum(dat.time(end));
era5.time = ttime(ttime>t1 & ttime<t2);
era5.slp = nan*era5.time;
era5.rh = nan*era5.time;
for i = 1:length(era5.time)
    time_idx = ttime == era5.time(i);
    [~,idx2] = min(abs(dat.time - era5.time(i)),[],'omitnan');
    [~,lon_idx] = min(abs(lon-dat.longitude(idx2)),[],'omitnan');
    [~,lat_idx] = min(abs(lat-dat.latitude(idx2)),[],'omitnan');
    era5.slp(i) = slp(lon_idx,lat_idx,time_idx); 
    era5.rh(i) = RH(lon_idx,lat_idx,time_idx);
    era5.ta(i) = t2mC(lon_idx,lat_idx,time_idx);
end
save([path_name,'sunfish_era5_data.mat'],'era5')

% extract pearldiver data
path_name = './../mat_files/'
var_name  = 'pearldiver_data'
load(fullfile(path_name,[var_name,'_clean.mat']))
dat = pearldiver;  clear pearldiver
dat.time = datenum(dat.time);

t1 = datenum(dat.time(1));
t2 = datenum(dat.time(end));
era5.time = ttime(ttime>t1 & ttime<t2);
era5.slp = nan*era5.time;
era5.rh = nan*era5.time;
for i = 1:length(era5.time)
    time_idx = ttime == era5.time(i);
    [~,idx2] = min(abs(dat.time - era5.time(i)),[],'omitnan');
    [~,lon_idx] = min(abs(lon-dat.longitude(idx2)));
    [~,lat_idx] = min(abs(lat-dat.latitude(idx2)));
    era5.slp(i) = slp(lon_idx,lat_idx,time_idx); 
    era5.rh(i) = RH(lon_idx,lat_idx,time_idx);
    era5.ta(i) = t2mC(lon_idx,lat_idx,time_idx);
end
save([path_name,'pearldiver_era5_data.mat'],'era5')

