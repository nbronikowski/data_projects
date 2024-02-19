clear; close all; clc; clear all;
flist = dir('*.nc');
flist(2)=[];
ttime = double(ncread(flist.name,'time'))/24+datenum(1900,1,1);
lon  = double(ncread(flist.name,'longitude'));
lat  = double(ncread(flist.name,'latitude'));
v10  = double(ncread(flist.name,'v10')); % mean sea level pressure Pa
u10  = double(ncread(flist.name,'u10'));  % surface pressure Pa
sst  = double(ncread(flist.name,'sst'))+273.15;
[long,latg]=meshgrid(lon,lat);


%% extract sunfish data
path_name = './../mat_files/'
var_name  = 'sunfish_data';
load(fullfile(path_name,[var_name,'_oxy_qc.mat']));
dat = sunfish; clear sunfish
dat.time = datenum(dat.time);


t1 = datenum(dat.time(1));
t2 = datenum(dat.time(end));
era5.time = ttime(ttime>t1 & ttime<t2);
era5.ekd = nan*era5.time;
era5.wd  = nan*era5.time;
era5.u10 = nan*era5.time;
era5.v10 = nan*era5.time;
era5.sst = nan*era5.time;

for i = 1:length(era5.time)
    time_idx = ttime == era5.time(i);
    [~,idx2] = min(abs(dat.time - era5.time(i)),[],'omitnan');
    [~,lon_idx] = min(abs(lon-dat.longitude(idx2)));
    [~,lat_idx] = min(abs(lat-dat.latitude(idx2)));
    
    [~,~,temp_Wd,temp_EKd] = ekman(latg,long,u10(:,:,time_idx)',v10(:,:,time_idx)','rho',1027.65);

    era5.wd(i)  = temp_Wd(lat_idx,lon_idx,1);
    era5.ekd(i) = temp_EKd(lat_idx,lon_idx,1);
    era5.u10(i) = u10(lon_idx,lat_idx,time_idx); 
    era5.v10(i) = v10(lon_idx,lat_idx,time_idx);
    era5.sst(i) = sst(lon_idx,lat_idx,time_idx);
end

save([path_name,'sunfish_era5_wind_ekman_data.mat'],'era5')



