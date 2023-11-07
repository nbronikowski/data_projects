clear; clc

% Load model CMEMS
load ADCP_inversion_result.mat

% Define the time period, location, and depth information
startDate = '2023-06-01'; % Adjust as necessary
endDate = '2023-06-13'; % Adjust as necessary

set(gca,'xlim',[-54.9 -50.5],'ylim',[46.7 50.9]);

lat_min = 46.5;
lat_max = 51.0;
lon_min = -55.0;
lon_max = -50.0;
depth_min = 0.493; % in pressure levels, not meters
depth_max = 350.00; % in pressure levels, not meters

% OPeNDAP URL
url = 'https://nrt.cmems-du.eu/thredds/dodsC/cmems_mod_glo_phy-cur_anfc_0.083deg_PT6H-i';

% Find the indices corresponding to the desired latitude and longitude ranges
info = ncinfo(url);
lat = ncread(url, 'latitude');
lon = ncread(url, 'longitude');
depth = ncread(url, 'depth');
time = ncread(url, 'time');

latIndices = find(lat >= lat_min & lat <= lat_max);
lonIndices = find(lon >= lon_min & lon <= lon_max);
depthIndices = find(depth >= depth_min & depth <= depth_max);

% Convert dates to the format "hours since 1950-01-01 00:00:00"
referenceDate = datenum('1950-01-01 00:00:00', 'yyyy-mm-dd HH:MM:SS');
time = time/24+referenceDate;

timeStart = datenum(startDate, 'yyyy-mm-dd'); % Convert start date to time index
timeEnd = datenum(endDate, 'yyyy-mm-dd') ; % Convert end date to time index

% Find the indices corresponding to the desired time range
timeIndices = find(time >= timeStart & time <= timeEnd);
start = [1+lonIndices(1), 1+latIndices(1), 1+min(depthIndices), 1+timeIndices(1)]; % [longitude, latitude, depth, time]
count = [length(lonIndices), length(latIndices), length(depthIndices), length(timeIndices)]; % [longitude, latitude, depth, time]

% count = [lonIndices(end), latIndices(end), timeIndices(end), depthIndices(end)];
% Read the data subset
u_data = ncread(url, 'uo', start, count);
v_data = ncread(url, 'vo', start, count);

m_time= time(timeIndices);
m_lon = lon(lonIndices);
m_lat = lat(latIndices);
m_depth= depth(depthIndices);


% Display the size of the retrieved data subset
gg_lon = out.glider_lon;
gg_lat = out.glider_lat;
gg_time= scale_var(out.seg_mid_time,1/4); % make it 6 hr timeseries
[g_time,idx] = unique(gg_time);
g_lon = gg_lon(idx);
g_lat = gg_lat(idx);
idx = isnan(g_time);
g_lon(idx)=[];
g_lat(idx)=[];
g_time(idx)=[];

u_model_glider = NaN(length(m_depth),length(g_time));
v_model_glider = NaN(length(m_depth),length(g_time));

for i = 1:length(g_time)
    [~,idTME] = nanmin(abs(m_time-g_time(i)));
    [~,idLON] = nanmin(abs(m_lon-g_lon(i)));
    [~,idLAT] = nanmin(abs(m_lat-g_lat(i)));
    u_model_glider(:,i) = squeeze(u_data(idLON,idLAT,:,idTME));
    v_model_glider(:,i) = squeeze(v_data(idLON,idLAT,:,idTME));
end

out.u_model = u_model_glider;
out.v_model = v_model_glider;
out.time_model = g_time;
out.lon_model = g_lon;
out.lat_model = g_lat;
out.depth_model = m_depth;


%% PLOT ADCP and CMEMS Model
figure()

subplot(221)
pcolor(out.seg_mid_time,-out.z_grid,out.ad2cp_u_z); shading interp
colormap(gca,cmocean('balance'))
caxis([-0.5 0.5])
ylim([-out.z_grid(end) out.z_grid(1)])
xlim([nanmin(out.seg_mid_time) nanmax(out.seg_mid_time)])
datetick('x','dd.mm','keeplimits')
title('East-West Currents')
ylabel('Glider ADCP Depth (m)')
formatplot

subplot(222)
pcolor(out.seg_mid_time,-out.z_grid,out.ad2cp_u_z); shading interp
colormap(gca,cmocean('balance'))
caxis([-0.5 0.5])
ylim([-out.z_grid(end) out.z_grid(1)])
xlim([nanmin(out.seg_mid_time) nanmax(out.seg_mid_time)])
datetick('x','dd.mm','keeplimits')
title('North-South Currents')
formatplot

subplot(223)
pcolor(out.time_model,-out.depth_model,out.u_model); shading interp
colormap(gca,cmocean('balance'))
caxis([-0.5 0.5])
ylim([-out.z_grid(end) out.z_grid(1)])
xlim([nanmin(out.seg_mid_time) nanmax(out.seg_mid_time)])
datetick('x','dd.mm','keeplimits')
ylabel('Model Depth (m)')
formatplot

subplot(224)
pcolor(out.time_model,-out.depth_model,out.v_model); shading interp
colormap(gca,cmocean('balance'))
caxis([-0.5 0.5])
ylim([-out.z_grid(end) out.z_grid(1)])
xlim([nanmin(out.seg_mid_time) nanmax(out.seg_mid_time)])
datetick('x','dd.mm','keeplimits')
formatplot

save_figure(gcf,'./plots/ADCP_CMEMS',[7.5 6],'.png','300')
