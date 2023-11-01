clear; close all; clc

fname = 'unit_334_delayed_trajectory.nc';
fpath = '';
fname = [fpath,fname];
gdat.time = ncread(fname,'time');

% Science Sensors
gdat.temp = ncread(fname,'temperature');
gdat.salt = ncread(fname,'salinity');
gdat.pres = ncread(fname,'pressure');
gdat.cond = ncread(fname,'conductivity');
gdat.depth = ncread(fname,'depth');
%gdat.chla = ncread(fname,'chlorophyll_a');
gdat.bb700_units = ncread(fname,'/glider_record/sci_flbb_bb_units');
gdat.bb700_sig = ncread(fname,'/glider_record/sci_flbb_bb_sig');
gdat.bb700_ref = ncread(fname,'/glider_record/sci_flbb_bb_ref');
gdat.chlor_units = ncread(fname,'/glider_record/sci_flbb_chlor_units');
gdat.chlor_sig = ncread(fname,'/glider_record/sci_flbb_chlor_sig');
gdat.chlor_ref = ncread(fname,'/glider_record/sci_flbb_chlor_ref');

% Navigation
gdat.lon  = ncread(fname,'lon_qc');
gdat.lat  = ncread(fname,'lat_qc');
gdat.prof_time = datenum(datetime(ncread(fname,'profile_time'),'ConvertFrom','posixtime'));
gdat.prof_idx = ncread(fname,'profile_index');
gdat.prof_dir = ncread(fname,'profile_direction');
gdat.m_depth    = ncread(fname,'glider_record/m_depth');
gdat.m_depth_state = ncread(fname,'glider_record/m_depth_state');
gdat.m_water_vx = ncread(fname,'glider_record/m_water_vx');
gdat.m_water_vy = ncread(fname,'glider_record/m_water_vy');
gdat.m_final_water_vx = ncread(fname,'glider_record/m_final_water_vx');
gdat.m_final_water_vy = ncread(fname,'glider_record/m_final_water_vy');
gdat.m_initial_water_vx = ncread(fname,'glider_record/m_initial_water_vx');
gdat.m_initial_water_vy = ncread(fname,'glider_record/m_initial_water_vy');
gdat.x_dr_state = ncread(fname,'glider_record/x_dr_state');
gdat.u_dac = ncread(fname,'u');
gdat.v_dac = ncread(fname,'v');

idnan = ~isnan(gdat.lon) | ~isnan(gdat.lat) ;
gdat.lon=interp1gap(gdat.time(idnan),gdat.lon(idnan),gdat.time,3600);
gdat.lat=interp1gap(gdat.time(idnan),gdat.lat(idnan),gdat.time,3600);

idnan = ~isnan(gdat.x_dr_state);
gdat.x_dr_state=interp1(gdat.time(idnan),gdat.x_dr_state(idnan),gdat.time,'nearest');

gdat.m_heading = ncread(fname,'glider_record/m_heading');
gdat.m_roll = ncread(fname,'glider_record/m_roll');
gdat.m_pitch = ncread(fname,'glider_record/m_pitch');

gdat.water_depth = ncread(fname,'glider_record/m_water_depth');

%% GET RID OF DUPLICATE TIMES DUE TO SMALL SHIFTS IN TIME
gdat.timeDateNum = datenum(datetime(gdat.time,'ConvertFrom','posixtime'));
[utimes,uidx] = unique(gdat.timeDateNum);
fn = fieldnames(gdat);
for i = 1:length(fn)
    temp = gdat.(fn{i});
    gdat.(fn{i}) = temp(uidx); %#ok<*FNDSB> 
end

% select only data from deployment
t1 = datenum(2023,08,03);
t2 = datenum(2023,09,24);
fn = fieldnames(gdat);
idx = find(gdat.timeDateNum>=t1 & gdat.timeDateNum<=t2);
for i = 1:length(fn)
    gdat.(fn{i}) = gdat.(fn{i})(idx); %#ok<*FNDSB> 
end

save([fname,'.mat'],'gdat')

