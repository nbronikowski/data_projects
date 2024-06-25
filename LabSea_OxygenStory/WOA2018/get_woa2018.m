clear; close all; clc

flist=dir('*nc')

for i = 1:length(flist)
    if i == 1
        woa2018.depth=ncread(flist(i).name,'depth');
        woa2018.lon=ncread(flist(i).name,'lon');
        woa2018.lon(woa2018.lon>180) = woa2018.lon(woa2018.lon > 180)-360;
        woa2018.lat=ncread(flist(i).name,'lat');
        depth_idx = woa2018.depth<20;
    end
    woa2018.time(i)=datenum(2018,i,1);
    temp_obj_oxy_sat=ncread(flist(i).name,'O_an');
    temp_oxy_sat=ncread(flist(i).name,'O_mn');
    temp_obj_oxy_sat(abs(temp_obj_oxy_sat)>1e5)=NaN;
    temp_oxy_sat(abs(temp_oxy_sat)>1e5)=NaN;
    woa2018.obj_oxy_sat(:,:,i) = squeeze(mean(temp_obj_oxy_sat(:,:,...
        depth_idx),3,'omitnan'))';
    woa2018.mn_oxy_sat(:,:,i) = squeeze(mean(temp_oxy_sat(:,:,...
        depth_idx),3,'omitnan'))'; % transform to lat lon representation
end

save('woa_oxy_sat_2018.mat','woa2018');

pcolor(woa2018.lon,woa2018.lat,woa2018.obj_oxy_sat(:,:,1))

