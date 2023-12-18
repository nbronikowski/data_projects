clear; clc; close all;

path_name = './mat_files/'
var_name  = 'sunfish_data'
load(fullfile(path_name,[var_name,'_clean.mat']))
dat = sunfish; 

% load ERA5 hourly single level data
load('./mat_files/sunfish_era5_data.mat')

slp = interp1(era5.time,era5.slp,dat.gridded.timeg);
rh  = interp1(era5.time,era5.rh, dat.gridded.timeg);
ta  = interp1(era5.time,era5.ta, dat.gridded.timeg);
pH2O  = vpress(nanmean(dat.gridded.salinity(1:5,:),1),nanmean(dat.gridded.temperature(1:5,:),1)); % atm top 5 m

% Compute Partial Pressure Oxygen Air and Wet
pO2_opt = O2ctoO2p(dat.gridded.oxygen_raw,dat.gridded.temperature,...
    dat.gridded.salinity,dat.gridded.pressure); % mbar Optode derived

pO2_air = (slp/101325-pH2O.*rh)*0.20946*1013.25; % mbar ERA5 derived


psat = satvap(ta,slp*0.01);
pO2_ref = 0.20946*(1013.25-psat);               % mbar Theoretical in-air

%% Next Extract Surface Measuring Intervals

pO2wet_1m = nanmean(pO2_opt(1:2,:),1); % surface O2
pO2wet_2m = nanmean(pO2_opt(3:5,:),1); % subsurface
depth_1m  = nanmean(dat.gridded.pressure(1:2,:),1);


% This is getting under/super saturation
d_O2wet_a = (pO2wet_1m ./pO2_ref)-1; % glider near  air 1m
d_O2air_a = (pO2_air ./pO2_ref)-1; % era5 data 
d_O2wet_w2m = (pO2wet_2m ./ pO2_ref)-1;

%% Carry Over Effect:
% From github ARGO Canada DMQC
x1 = pO2wet_1m - pO2_air;
y1 = pO2wet_2m - pO2_air;
x1 = x1(:); y1 = y1(:);

id = ~isnan(x1) & ~isnan(y1);
x1 = x1(id);
y1 = y1(id);
X = [x1, ones(length(x1), 1)];  % Adding a column of ones for intercept
c = X \ y1;
c = c(1);

O2_gains = ((1-c)*pO2_air)./(pO2wet_1m - c*pO2wet_2m);

% plot(dat.gridded.timeg,1./O2_gains,'.')

% This is how I would estimate the gain ref vs measured to adjust
% O2_gains= (pO2air_1m./ pO2wet_1m); % simple absolute gain based on pO2
% air should be pO2 glider at surface
O2_gains_o= (d_O2air_a+1)./(d_O2wet_a+1); 

%% Stuff from Nicholson 2017 -- trying to recreate plots
% d_O2wet_a = (pO2wet_1m ./dat.pO2ref)-1; % glider near  air 1m
% d_O2air_a = (pO2air_1m ./dat.pO2ref)-1; % era5 data 
% d_O2wet_w2m = (pO2wet_2m ./dat.pO2ref)-1;
% d_O2wet_w100m = (pO2wet_100m./dat.pO2ref)-1;

% idnan = isnan(d_O2wet_a) | isnan(d_O2wet_w100m) | isnan(d_O2wet_w2m);
% d_O2wet_w100m(idnan)=[];
% d_O2wet_a(idnan)=[];
% d_O2wet_w2m(idnan)=[];

tthresh = datenum(2022,04,01);
id = find(dat.gridded.timeg(~isnan(dat.gridded.timeg)) < tthresh);
id2 = find(dat.gridded.timeg(~isnan(dat.gridded.timeg)) > tthresh);


% Create a figure
figure;

% Create a tiled layout with 1 row and 2 columns
t = tiledlayout(2, 2);

% Create the first tile
nexttile([1 2]);
hold on;
h1 = plot(dat.gridded.timeg, O2_gains, '.');
h2 = plot(dat.gridded.timeg, dat.gridded.timeg*0 + nanmedian(O2_gains(id)), '-b');
h3 = plot(dat.gridded.timeg, dat.gridded.timeg*0 + nanmean(O2_gains(id)), '-r');
h4 = plot(dat.gridded.timeg, dat.gridded.timeg*0 + nanmean(O2_gains(id)) - std(O2_gains(id), [], "all", "omitnan"), '--k');
plot(dat.gridded.timeg, dat.gridded.timeg*0 + nanmean(O2_gains(id)) + std(O2_gains(id), [], "all", "omitnan"), '--k');
plot([tthresh tthresh], [0.9 1.2], ':', 'Color', [.6 .6 .6]);
legend([h1 h2 h3 h4], {'gains', ...
    ['median (',num2str(round(nanmedian(O2_gains(id)),3)),')'], ...
    ['mean (',num2str(round(nanmean(O2_gains(id)),3)),')'], ...
    ['\pm 1\sigma (',num2str(round(nanstd(O2_gains(id)),3)),')'],...
    }, 'Location', 'best');
ylabel('gain (\Delta O^{air}_{2,a} / \Delta O^{wet}_{2,a})');
xlim([datenum(2022, 01, 01) datenum(2022, 05, 30)]);
datetick('x', 'dd-mmm-yy');
title('(a)','HorizontalAlignment','left')
formatplot;

% Create the second tile
nexttile;
hold on;
h = histogram(O2_gains(id));
h2 = histogram(O2_gains(id2), 'FaceColor', 'm', 'BinWidth', h.BinWidth);
h3 = plot([nanmean(O2_gains(id)) nanmean(O2_gains(id))], [0 100], '-r', 'LineWidth', 2);
h4 = plot([nanmedian(O2_gains(id)) nanmedian(O2_gains(id))], [0 100], '-b', 'LineWidth', 2);
xlabel('gain factor \Delta O^{air}_{2,a} / \Delta O^{wet}_{2,a}');
ylabel('count');
title('(b)','HorizontalAlignment','left')
legend('gains (t<01-04-2022)', 'gains (t>01-04-2022)', 'mean', 'median', 'Location', 'best');
formatplot;


% Format date tick labels
nexttile; hold on
plot(d_O2wet_w2m*100+100,d_O2wet_a*100+100,'*');
plot([82 95],[82 95],':k')
xlim([82 95])
ylim([82 95])
ylabel('O^{wet,a}_2 / %')
xlabel('O^{wet, 2m}_2 / %')
formatplot
title('(c)','HorizontalAlignment','left')
legend('glider data','1:1 fit','location','best')


t.TileSpacing='compact';
t.Padding='compact';
subtitle(t,'Sunfish Glider and ERA5 derived Surface pO_2 Gains (Carry-Over Effect)')
save_figure(gcf,['./plots/sunfish_era5_optode_gains_carry_over'],[7.5 5],['.png'],'300')

%% Adjust Raw Oxygen using the median gain before May 1 calculated from the ERA5 pO2 comparison
tthresh = datenum(2022,04,01);
id = find(dat.gridded.timeg < tthresh);
sunfish.adjusted_oxygen_concentration = dat.raw_oxygen_concentration*nanmean(O2_gains(id));
sunfish.gridded.oxygen_adjusted = dat.gridded.oxygen_raw*nanmean(O2_gains(id));

% save result
save(fullfile(path_name,[var_name,'_oxy_qc.mat']),'sunfish','-v7.3')

