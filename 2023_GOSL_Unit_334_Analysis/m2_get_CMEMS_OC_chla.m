clear; clc

% Define the time period, location, and depth information
startDate = '2023-08-10'; % Adjust as necessary
endDate = '2023-09-02'; % Adjust as necessary
lat_min = 49;
lat_max = 51.5;
lon_min = -60;
lon_max = -57.5;

% OPeNDAP URL
url = 'https://my.cmems-du.eu/thredds/dodsC/cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D';

% Find the indices corresponding to the desired latitude and longitude ranges
info = ncinfo(url);
lat = double(ncread(url, 'lat'));
lon = double(ncread(url, 'lon'));
time = double(ncread(url, 'time'));

latIndices = find(lat >= lat_min & lat <= lat_max);
lonIndices = find(lon >= lon_min & lon <= lon_max);

% Convert dates to the format days "since 1900-01-01"
referenceDate = datenum(1900,1,1)-1;
time = time+referenceDate;

timeStart = datenum(startDate, 'yyyy-mm-dd'); % Convert start date to time index
timeEnd = datenum(endDate, 'yyyy-mm-dd') ; % Convert end date to time index

% Find the indices corresponding to the desired time range
timeIndices = find(time >= timeStart & time <= timeEnd);
start = [1+lonIndices(1), 1+latIndices(1), 1+timeIndices(1)]; % [longitude, latitude, depth, time]
count = [length(lonIndices), length(latIndices), length(timeIndices)]; % [longitude, latitude, depth, time]

% count = [lonIndices(end), latIndices(end), timeIndices(end), depthIndices(end)];
% Read the data subset
chla = ncread(url, 'CHL', start, count);
chla_err = ncread(url, 'CHL_uncertainty', start, count);


cm.time  = time(timeIndices);
cm.lon  = lon(lonIndices);
cm.lat  = lat(latIndices);
cm.chla = chla;
cm.chla_err = chla_err;

save('CMEMS_multi_chla.mat','cm');

%% Plot
clear; close all; clc

load CMEMS_multi_chla.mat

PanFac = 4; % we want 4 columns

% Create a tiled layout with desired dimensions
t = tiledlayout(PanFac, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

k = 1;

for i = 11:length(cm.time)-2
    
    % Use nexttile instead of subplot
    ax = nexttile; hold on
    imagescn(cm.lon,cm.lat,cm.chla(:,:,i)')
    borders('Canada','facecolor',[.6 .6 .6])
    colormap(ax,cmocean('algae'))
    caxis([0 2]);
    title(datestr(cm.time(i),'yyyy-mmm-dd'))
    ylim([lat_min lat_max])
    xlim([lon_min lon_max])
    if i == 22
        cb = colorbar;
        cb.Layout.Tile = 'east'; % To position the colorbar on the east side of the tiled layout
    end
end

