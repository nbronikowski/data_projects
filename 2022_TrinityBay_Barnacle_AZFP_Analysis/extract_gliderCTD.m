clear; close all; clc

fpath = './glider_data/';
fname = 'barnacle_trinity_bay_mission_L1.nc';

var_list = {...
'time'
'time_unix'
'time_ctd'
'lat'
'lon'
'wpt_lat'
'wpt_lon'
'distance_over_ground'
'pitch'
'roll'
'heading'
'water_u'
'water_v'
'profile_index'
'profile_direction'
'depth'
'temperature'
'pressure'
'conductivity'
'salinity'
'abs_salinity'
'density'};

dat = [];
for i = 1:length(var_list)
    dat.(var_list{i})=ncread([fpath,fname],var_list{i});
end

uprofs = unique(dat.profile_index(~isnan(dat.profile_index)));
uprofs = uprofs(mod(uprofs,1)==0);

profile_index = dat.profile_index;
profile_direction = dat.profile_direction;

profile_time  = [];
dcast = struct();
for j =1:length(var_list)
    dcast.(var_list{j})=[];
end
% We want profile time and we only want downcasts
for i = 1:length(uprofs)
    idx = find(profile_index==uprofs(i) & profile_direction==1);
    profile_time = [profile_time;mean(dat.time(idx),'omitnan')*ones(size(dat.time(idx)))];
    for j = 1:length(var_list)
        dcast.(var_list{j}) = [dcast.(var_list{j});dat.(var_list{j})(idx)];
    end
end
dcast.profile_time = profile_time;


dcast.sound_speed = sw_svel(dcast.salinity,dcast.temperature,dcast.pressure);
dat.sound_speed = sw_svel(dat.salinity,dat.temperature,dat.pressure);
save('glider_ctd_data.mat','dat')
save('glider_ctd_dcast.mat','dcast');