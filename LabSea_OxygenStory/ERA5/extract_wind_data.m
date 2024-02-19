clear; clc; close all

flist = dir('*.nc');
flist(2)=[];
ttime = double(ncread(flist.name,'time'))/24+datenum(1900,1,1);
lon  = double(ncread(flist.name,'longitude'));
lat  = double(ncread(flist.name,'latitude'));
v10  = double(ncread(flist.name,'v10')); % mean sea level pressure Pa
u10  = double(ncread(flist.name,'u10'));  % surface pressure Pa
sst  = double(ncread(flist.name,'sst'))+273.15;

% extract sunfish data
path_name = './../mat_files/'
var_name  = 'sunfish_data'
load(fullfile(path_name,[var_name,'_clean.mat']))
dat = sunfish;  clear sunfish
dat.time = datenum(dat.time);

t1 = datenum(dat.time(1));
t2 = datenum(dat.time(end));
era5.time = ttime(ttime>t1 & ttime<t2);
era5.u10 = nan*era5.time;
era5.v10 = nan*era5.time;
era5.sst = nan*era5.time;
for i = 1:length(era5.time)
    time_idx = ttime == era5.time(i);
    [~,idx2] = min(abs(dat.time - era5.time(i)),[],'omitnan');
    [~,lon_idx] = min(abs(lon-dat.longitude(idx2)));
    [~,lat_idx] = min(abs(lat-dat.latitude(idx2)));
    era5.u10(i) = u10(lon_idx,lat_idx,time_idx); 
    era5.v10(i) = v10(lon_idx,lat_idx,time_idx);
    era5.sst(i) = sst(lon_idx,lat_idx,time_idx);
end
save([path_name,'sunfish_era5_wind_data.mat'],'era5')

% extract pearldiver data
path_name = './../mat_files/'
var_name  = 'pearldiver_data'
load(fullfile(path_name,[var_name,'_clean.mat']))
dat = pearldiver;  clear pearldiver
dat.time = datenum(dat.time);

t1 = datenum(dat.time(1));
t2 = datenum(dat.time(end));
era5.time = ttime(ttime>t1 & ttime<t2);
era5.u10 = nan*era5.time;
era5.v10 = nan*era5.time;
era5.sst = nan*era5.time;
for i = 1:length(era5.time)
    time_idx = ttime == era5.time(i);
    [~,idx2] = min(abs(dat.time - era5.time(i)),[],'omitnan');
    [~,lon_idx] = min(abs(lon-dat.longitude(idx2)));
    [~,lat_idx] = min(abs(lat-dat.latitude(idx2)));
    era5.u10(i) = u10(lon_idx,lat_idx,time_idx); 
    era5.v10(i) = v10(lon_idx,lat_idx,time_idx);
    era5.sst(i) = sst(lon_idx,lat_idx,time_idx);
end
save([path_name,'pearldiver_wind_era5_data.mat'],'era5')

