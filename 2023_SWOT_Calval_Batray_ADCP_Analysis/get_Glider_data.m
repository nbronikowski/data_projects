clear; close all; clc

addpath('./glider_data/')
fname = 'batray_2023-06-01_delayed_trajectory.nc';
% fname = [fpath,fname];

% Navigation Data
gdat.time = ncread(fname,'time');
gdat.lon  = ncread(fname,'lon');
gdat.lat  = ncread(fname,'lat');

gdat.m_heading = ncread(fname,'glider_record/m_heading');
gdat.m_roll = ncread(fname,'glider_record/m_roll');
gdat.m_pitch = ncread(fname,'glider_record/m_pitch');

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


% Interpolation U,V, LON, LAT
gdat.u_dac = ncread(fname,'u');
gdat.v_dac = ncread(fname,'v');
idnan = ~isnan(gdat.u_dac) | ~isnan(gdat.v_dac) ;
gdat.u_dac=interp1gap(gdat.time(idnan),gdat.u_dac(idnan),gdat.time,3600);
gdat.v_dac=interp1gap(gdat.time(idnan),gdat.v_dac(idnan),gdat.time,3600);

idnan = ~isnan(gdat.lon) | ~isnan(gdat.lat) ;
gdat.lon=interp1gap(gdat.time(idnan),gdat.lon(idnan),gdat.time,3600);
gdat.lat=interp1gap(gdat.time(idnan),gdat.lat(idnan),gdat.time,3600);

idnan = ~isnan(gdat.x_dr_state);
gdat.x_dr_state=interp1(gdat.time(idnan),gdat.x_dr_state(idnan),gdat.time,'nearest');


% water depth
gdat.water_depth = ncread(fname,'glider_record/m_water_depth');
gdat.water_depth(gdat.water_depth<20)=NaN;
idnan = find(~isnan(gdat.water_depth));
temp = gdat.water_depth(idnan); idPeaks = abs([0;diff(temp)])>5; temp(idPeaks)=NaN;
gdat.water_depth(idnan)=temp;
idnan = find(~isnan(gdat.water_depth));
gdat.water_depth=interp1gap(gdat.time(idnan),gdat.water_depth(idnan),gdat.time,3600);

% CTD
gdat.temp = ncread(fname,'temperature');
gdat.ctemp= ncread(fname,'conservative_temperature');
gdat.salt = ncread(fname,'salinity');
gdat.cond = ncread(fname,'conductivity')/10; % forgot to fix this in glider processing.
gdat.saltSA=ncread(fname,'absolute_salinity');
gdat.pres = ncread(fname,'pressure');
gdat.depth = ncread(fname,'depth');
gdat.c_stp= sw_svel(gdat.salt,gdat.temp,gdat.pres);

% Oxygen
gdat.oxygen_concentration = ncread(fname,'oxygen_concentration');
gdat.oxygen_calphase = ncread(fname,'glider_record/sci_oxy4_calphase');
gdat.oxygen_temp = ncread(fname,'oxygen_sensor_temperature');

% Oxygen Recalculation
foilcoef = [3.10365E-03 1.32185E-04 2.60868E-06 2.07958E02 -2.36738E-01 ...
    -4.25867E01 4.07417E00];
modeltype = 'uchidaAADI';
gdat.oxygen_phase_recalc = optcalcO2(gdat.temp,gdat.oxygen_calphase,foilcoef,...
    modeltype,gdat.salt,1013.25,gdat.pres);

%% GET RID OF DUPLICATE TIMES DUE TO SMALL SHIFTS IN TIME
gdat.timeDateNum = datenum(datetime(gdat.time,'ConvertFrom','posixtime'));
[utimes,uidx] = unique(gdat.timeDateNum);
fn = fieldnames(gdat);
for i = 1:length(fn)
    temp = gdat.(fn{i});
    gdat.(fn{i}) = temp(uidx); %#ok<*FNDSB> 
end

%% Remove data outside range
t1 = datenum(2023,06,01,0,0,0);
t2 = datenum(2023,06,27,0,0,0);
idx = find(gdat.timeDateNum>=t1 & gdat.timeDateNum<=t2);
for i = 1:length(fn)
    gdat.(fn{i}) = gdat.(fn{i})(idx); %#ok<*FNDSB> 
end

gdat.lon(gdat.lat>80)=NaN;
gdat.lat(gdat.lat>80)=NaN;

%% PLOT TRACK
figure()
plot(gdat.lon,gdat.lat,'.r','MarkerSize',10); hold on; 
borders('Canada','facecolor','g');
% quiver(gdat.lon,gdat.lat,gdat.m_initial_water_vx,gdat.m_initial_water_vy,150,'filled','Color','r');
% quiver(gdat.lon,gdat.lat,gdat.m_final_water_vx,gdat.m_final_water_vy,150,'filled','Color','b');
quiver(gdat.lon,gdat.lat,gdat.m_final_water_vx,gdat.m_final_water_vy,150,'filled','Color','m');
set(gca,'xlim',[-54.9 -50.5],'ylim',[46.7 50.9]);
xlim = get(gca,'XLim'); ylim = get(gca,'Ylim');
[x,y,z]=load_bathy(ylim,xlim);
title('Glider Deployment in Support of SWOT Calval Mission')
contour(x,y,z','color','k','levellist',[-1000,-300,-200,-100,-50],'showtext','on');
save_figure(gcf,['./plots/glider_track'],[7.5 5],'.png','300')


% process_CTD
% process_OPTODE

save('glider_data.mat','gdat')


