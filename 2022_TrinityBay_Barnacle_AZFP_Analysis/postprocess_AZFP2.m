clear; clc; close all;

load('azfp_out.mat');
load('glider_ctd_dcast.mat');

% Merge and Process AZFP - glider record
azfp_time = unique(azfp(1).Date,'stable') + 1e-6;
tgap = 300/86400;
azfp_depth = interp1gap(dcast.time, dcast.depth, azfp_time, tgap);
azfp_pidx = interp1(dcast.time, dcast.profile_index, azfp_time, 'nearest');
azfp_lat = interp1gap(dcast.time, dcast.lat, azfp_time, tgap);
azfp_lon = interp1gap(dcast.time, dcast.lon, azfp_time, tgap);
ind = isnan(azfp_depth);
azfp_depth(ind) = [];
azfp_pidx(ind) = [];
azfp_lat(ind) = [];
azfp_lon(ind) = [];
azfp_time(ind) = [];

%% Try to correct tilt using ADCP code using orientation for the 3 rd 
%  straight down case and average bins in same depth for every profile




% initialize the Channel 1 -> 4
echoSv = azfp(4).Sv;
echoSv(ind, :) = [];
bin_depth = azfp(3).Range(1, :);

% Find indices location of max counts for bottom detection
[~, loc] = max(echoSv(:, 10:end), [], 2);
ExtraY = 9; 
BottomLoc = loc + ExtraY;
[n, l] = max(BottomLoc);

% Initialize rotated matrices
rot_Sv = nan(size(echoSv, 1), n);
rot_depth = rot_Sv;
rot_time = rot_Sv;
tmp_depth = azfp_depth * ones(1, n) + ones(length(azfp_depth), 1) * bin_depth(1:n);
tmp_time = repmat(azfp_time, [1, n]);

% Rotate to have the floor signal at the bottom
for ii = 1:size(echoSv, 1)
    if BottomLoc(l) - BottomLoc(ii) + 1 < 1
        break;
    end
    range = BottomLoc(l) - BottomLoc(ii) + 1 : BottomLoc(l);
    rot_Sv(ii, range) = echoSv(ii, 1:BottomLoc(ii));
    rot_depth(ii, range) = tmp_depth(ii, 1:BottomLoc(ii));
    rot_time(ii, range) = tmp_time(ii, 1:BottomLoc(ii));
end

% Find the ends of the separate dives and assign each ping to their specific profile
pp = diff(datenum(azfp_time)) * 86400;
i_end = find(pp > 6); 
i_end = [i_end; length(azfp_time)]; % Add the last index
profile_nr = zeros(length(azfp_time), 1);

for i = 1:length(i_end)
    if i == 1
        start_idx = 1;
    else
        start_idx = i_end(i-1) + 1;
    end
    end_idx = i_end(i);
    profile_nr(start_idx:end_idx) = i;
end

% Combine data from each dive into a profile with .5 m bins
nr_dives = max(profile_nr);
dives_depth = 0:0.5:600;
dives_lat = nan(nr_dives, 1);
dives_lon = dives_lat;
dives_time = dives_lat;
dives_Sv = nan(nr_dives, length(dives_depth));


for i = 1:nr_dives
    ind_prof = find(profile_nr == i);
    dives_time(i) = nanmean(datenum(azfp_time(ind_prof)));
    dives_lat(i) = nanmean(azfp_lat(ind_prof));
    dives_lon(i) = nanmean(azfp_lon(ind_prof));
    tmp_Sv = rot_Sv(ind_prof, :);
    tmp_depth = rot_depth(ind_prof, :);
    
    for j = 1:length(dives_depth)
        ind_d = find(tmp_depth(:,:) > dives_depth(j) - .25 & tmp_depth(:,:) <= dives_depth(j) + .25);
        dives_Sv(i, j) = nanmean(tmp_Sv(ind_d));
    end
end


figure
set(gcf,'color','w');
h1=imagescn(dives_time,-dives_depth,dives_Sv')
datetick('x','mmm dd HH PM','keepticks');
xtickangle(45)
ylabel('Water depth (m)','FontSize',11)
xlabel('Time UTC (NST-3.5H UTC)','FontSize',11)
axis tight
cmocean('balance')
h2=colorbar;
ylabel(h2,'Volume backscatter Sv (dB re m^{-1})','FontSize',11);%,'Rotation',270.0)
set(gca,'YLim',[-600 0])
caxis([-90 -50])
