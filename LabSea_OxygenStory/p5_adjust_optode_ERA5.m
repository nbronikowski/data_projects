clear; clc; close all;

path_name = './mat_files/'
var_name  = 'pearldiver_data'
load(fullfile(path_name,[var_name,'_clean.mat']))
dat = pearldiver; 

% load woa2018 data and compute gains ...
load('./mat_files/pearldiver_era5_data.mat')

dat.time = datenum(dat.time);
idt = dat.time>datenum(2020,01,15) & dat.time<datenum(2020,05,15);
gtime = dat.time(idt);

id = ~isnan(dat.depth);
depth = interp1gap(dat.time,dat.depth,gtime,30/86400);

slp = interp1(era5.time,era5.slp,gtime);
rh  = interp1(era5.time,era5.rh, gtime);
ta  = interp1(era5.time,era5.ta, gtime);
pH2O  = vpress(dat.salinity(idt),dat.temperature(idt)); % atm

% Compute Partial Pressure Oxygen Air and Wet
pO2air = (slp/101325-pH2O.*rh)*0.20946*1013.25; % mbar
pO2wet = O2ctoO2p(dat.raw_oxygen_concentration(idt),dat.temperature(idt),...
    dat.salinity(idt),dat.pressure(idt));

psat = satvap(ta,slp*0.01);
pO2ref = 0.20946*(1013.25-psat); % mbar 

%% Next Extract Surface Measuring Intervals
ids=  depth<1; % picking surface points
pO2wet_1m = interp1gap(gtime(ids),pO2wet(ids),gtime,30/86400);
pO2air_1m = interp1gap(gtime(ids),pO2air(ids),gtime,30/86400);
depth_1m  = interp1gap(gtime(ids),depth(ids),gtime,30/86400);


uprofs = unique(dat.profile_index(~isnan(dat.profile_index)));
pO2wet_2m = [];
pO2wet_100m = [];
time_2m = [];
% extract ascending to surface pO2 measurements inside 1 - 2.5 m band
for i = 1:length(uprofs)
    prof_id = dat.profile_index(idt) == uprofs(i) ...
        & dat.profile_direction(idt) == 1;
    temp_pO2wet = pO2wet(prof_id);
    temp_depth  = depth(prof_id);
    temp_time   = gtime(prof_id);
    
    idz = temp_depth>1 & temp_depth<2.5;
    idz2 = temp_depth>100 & temp_depth<102;
    pO2wet_2m(i) = nanmean(temp_pO2wet(idz));
    pO2wet_100m(i) = nanmean(temp_pO2wet(idz2));
    time_2m(i) = nanmean(temp_time(idz));
    time_100m(i) = nanmean(temp_time(idz2));
end
id = ~isnan(pO2wet_2m);
pO2wet_2m = interp1(time_2m(id),pO2wet_2m(id),gtime);

id = ~isnan(pO2wet_100m);
pO2wet_100m = interp1(time_100m(id),pO2wet_100m(id),gtime);


% This is getting under/super saturation
d_O2wet_a = (pO2wet_1m ./pO2ref)-1; % glider near  air 1m
d_O2air_a = (pO2air_1m ./pO2ref)-1; % era5 data 
d_O2wet_w2m = (pO2wet_2m ./pO2ref)-1;

%% Carry Over Effect:
% From github ARGO Canada DMQC
x1 = pO2wet_1m - pO2air_1m;
y1 = pO2wet_2m - pO2air_1m;

id = ~isnan(x1) & ~isnan(y1);
x1 = x1(id);
y1 = y1(id);
X = [x1, ones(length(x1), 1)];  % Adding a column of ones for intercept
c = X \ y1;
c = c(1);

O2_gains = ((1-c)*pO2air_1m(id))./(pO2wet_1m(id) - c*pO2wet_2m(id));

plot(gtime(id),movmedian(O2_gains,30),'.')

% This is how I would estimate the gain ref vs measured to adjust
% O2_gains= (pO2air_1m./ pO2wet_1m); % simple absolute gain based on pO2
% air should be pO2 glider at surface
O2_gains_o= (d_O2air_a(id)+1)./(d_O2wet_a(id)+1); 

%% Stuff from Nicholson 2017 -- trying to recreate plots
d_O2wet_a = (pO2wet_1m ./dat.pO2ref)-1; % glider near  air 1m
d_O2air_a = (pO2air_1m ./dat.pO2ref)-1; % era5 data 
d_O2wet_w2m = (pO2wet_2m ./dat.pO2ref)-1;
d_O2wet_w100m = (pO2wet_100m./dat.pO2ref)-1;

idnan = isnan(d_O2wet_a) | isnan(d_O2wet_w100m) | isnan(d_O2wet_w2m);
d_O2wet_w100m(idnan)=[];
d_O2wet_a(idnan)=[];
d_O2wet_w2m(idnan)=[];

tthresh = datenum(2020,04,01);
id = find(dat.time(~isnan(dat.time)) < tthresh);
id2 = find(dat.time(~isnan(dat.time)) > tthresh);


% Create a figure
figure;

% Create a tiled layout with 1 row and 2 columns
t = tiledlayout(2, 2);

% Create the first tile
nexttile([1 2]);
hold on;
h1 = plot(dat.time(id), O2_gains, '.');
h2 = plot(dat.time, dat.time*0 + nanmedian(O2_gains(id)), '-b');
h3 = plot(dat.time, dat.time*0 + nanmean(O2_gains(id)), '-r');
h4 = plot(dat.time, dat.time*0 + nanmean(O2_gains(id)) - std(O2_gains(id), [], "all", "omitnan"), '--k');
plot(dat.time, dat.time*0 + nanmean(O2_gains(id)) + std(O2_gains(id), [], "all", "omitnan"), '--k');
plot([tthresh tthresh], [0.85 1.1], ':', 'Color', [.6 .6 .6]);
legend([h1 h2 h3 h4], {'gains', ...
    ['median (',num2str(round(nanmedian(O2_gains(id)),3)),')'], ...
    ['mean (',num2str(round(nanmean(O2_gains(id)),3)),')'], ...
    ['\pm 1\sigma (',num2str(round(nanstd(O2_gains(id)),3)),')'],...
    }, 'Location', 'best');
ylabel('gain (\Delta O^{air}_{2,a} / \Delta O^{wet}_{2,a})');
xlim([datenum(2020, 01, 01) datenum(2020, 05, 30)]);
datetick('x', 'dd-mmm-yy');
title('(a)','HorizontalAlignment','left')
formatplot;

% Create the second tile
nexttile;
hold on;
h = histogram(O2_gains(id));
h2 = histogram(O2_gains(id2), 'FaceColor', 'm', 'BinWidth', h.BinWidth);
h3 = plot([nanmean(O2_gains(id)) nanmean(O2_gains(id))], [0 2000], '-r', 'LineWidth', 2);
h4 = plot([nanmedian(O2_gains(id)) nanmedian(O2_gains(id))], [0 2000], '-b', 'LineWidth', 2);
xlabel('gain factor \Delta O^{air}_{2,a} / \Delta O^{wet}_{2,a}');
ylabel('count');
title('(b)','HorizontalAlignment','left')
legend('gains (t<01-04-2020)', 'gains (t>01-04-2020)', 'mean', 'median', 'Location', 'best');
formatplot;


% Format date tick labels
nexttile; hold on
% plot(d_O2wet_w2m*100+100,d_O2wet_a*100+100,'*');
plot(d_O2wet_w100m*100+100,d_O2wet_a*100+100,'^');
plot([90 105],[90 110],':k')
xlim([94 107])
ylim([94 107])
ylabel('O^{wet,a}_2 / %')
xlabel('O^{wet, 100m}_2 / %')
formatplot
title('(c)','HorizontalAlignment','left')
legend('glider data','1:1 fit','location','best')


t.TileSpacing='compact';
t.Padding='compact';
subtitle(t,'Pearldiver Glider and ERA5 derived Surface pO_2')
save_figure(gcf,['./plots/pearldiver_era5_optode_gains'],[7.5 5],['.png'],'300')

%% Adjust Raw Oxygen using the median gain before May 1 calculated from the ERA5 pO2 comparison
tthresh = datenum(2020,04,01);
id = find(dat.time(~isnan(dat.time)) < tthresh);

pearldiver.adjusted_oxygen_concentration = dat.raw_oxygen_concentration*nanmean(O2_gains(id));
pearldiver.gridded.oxygen_adjusted = dat.gridded.oxygen_raw*nanmean(O2_gains(id));

% save result
save(fullfile(path_name,[var_name,'_oxy_qc.mat']),'pearldiver','-v7.3')

