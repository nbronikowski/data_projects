clear; clc; close all;

path_name = './mat_files/'
var_name  = 'sunfish_data'
load(fullfile(path_name,[var_name,'_clean.mat']))
dat = sunfish;  clear sunfish

% load woa2018 data and compute gains ...
load('./mat_files/sunfish_era5_data.mat')

dat.time = datenum(dat.time);

id = ~isnan(dat.depth);
dat.depth = interp1gap(dat.time,dat.depth,dat.time,30/86400);

dat.slp = interp1(era5.time,era5.slp,dat.time);
dat.rh  = interp1(era5.time,era5.rh, dat.time);
dat.ta  = interp1(era5.time,era5.ta, dat.time);
dat.pH2O  = vpress(dat.salinity,dat.temperature); % atm

% Compute Partial Pressure Oxygen Air and Wet
dat.pO2air = (dat.slp/101325-dat.pH2O.*dat.rh)*0.20946*1013.25; % mbar
dat.pO2wet = O2ctoO2p(dat.raw_oxygen_concentration,dat.temperature,...
    dat.salinity,dat.pressure);

psat = satvap(dat.ta,dat.slp*0.01);
dat.pO2ref = 0.20946*(1013.25-psat); % mbar 

%% Next Extract Surface Measuring Intervals
ids=  dat.depth<1; % picking surface points
pO2wet_1m = interp1gap(dat.time(ids),dat.pO2wet(ids),dat.time,30/86400);
pO2air_1m = interp1gap(dat.time(ids),dat.pO2air(ids),dat.time,30/86400);
depth_1m  = interp1gap(dat.time(ids),dat.depth(ids),dat.time,30/86400);


uprofs = unique(dat.profile_index(~isnan(dat.profile_index)));
pO2wet_2m = [];
time_2m = [];
% extract ascending to surface pO2 measurements inside 1 - 2.5 m band
for i = 1:length(uprofs)
    prof_id = dat.profile_index == uprofs(i) ...
        & dat.profile_direction == -1;
    temp_pO2wet = dat.pO2wet(prof_id);
    temp_depth  = dat.depth(prof_id);
    temp_time   = dat.time(prof_id);
    
    idz = temp_depth>1 & temp_depth<2.5;
    pO2wet_2m(i) = nanmean(temp_pO2wet(idz));
    time_2m(i) = nanmean(temp_time(idz));
end
id = ~isnan(pO2wet_2m);
pO2wet_2m = interp1(time_2m(id),pO2wet_2m(id),dat.time);

% This is getting under/super saturation
d_O2wet_a = (pO2wet_1m ./dat.pO2ref)-1; % glider near  air 1m
d_O2air_a = (pO2air_1m ./dat.pO2ref)-1; % era5 data 

d_O2wet_w = (pO2wet_2m ./dat.pO2ref)-1;

% This is how I would estimate the gain ref vs measured to adjust
% O2_gains= (pO2air_1m./ pO2wet_1m); % simple absolute gain based on pO2
% air should be pO2 glider at surface

O2_gains= (d_O2air_a+1)./(d_O2wet_a+1); 

tthresh = datenum(2022,05,01);
id = find(dat.time(~isnan(dat.time)) < tthresh);


figure(); hold on
h1=plot(dat.time,O2_gains,'.'); 
h2=plot(dat.time,dat.time*0+nanmedian(O2_gains(id)),'-b');
h3=plot(dat.time,dat.time*0+nanmean(O2_gains(id)),'-r');
h4=plot(dat.time,dat.time*0+nanmean(O2_gains(id))-std(O2_gains(id),[],"all","omitnan"),'--k');
plot(dat.time,dat.time*0+nanmean(O2_gains(id))+std(O2_gains(id),[],"all","omitnan"),'--k');
% h6=polyplot(dat.time(id),O2_gains(id),'linestyle','-','color','b','linewidth',1.5);
plot([tthresh tthresh],[0.9 1.2],':','Color',[.6 .6 .6]);
legend([h1 h2 h3 h4],{'gains','median','mean','\pm 1\sigma'},'Location','best');
title('Sunfish Glider and ERA5 derived Surface pO_2')
ylabel('gain factor \Delta O^{air}_{2,a} / \Delta O^{wet}_{2,a} ' );
datetick('x','dd-mmm-yy')
formatplot
save_figure(gcf,['./plots/sunfish_era5_optode_gains'],[7.5 5],['.png'],'300')


%% Adjust Raw Oxygen using the median gain before May 1 calculated from the ERA5 pO2 comparison

dat.adjusted_oxygen_concentration = dat.raw_oxygen_concentration*nanmedian(O2_gains(id));
dat.gridded.oxygen_adjusted = dat.gridded.oxygen_raw*nanmedian(O2_gains(id));

% save result




%% Stuff from Nicholson 2017 -- trying to recreate plots 
id = isnan(d_O2wet_a) | isnan(d_O2wet_w);
d_O2wet_w(id)=[];
d_O2wet_a(id)=[];
d_O2air_a(id)=[];

figure()
subplot(121)
histogram(d_O2wet_a*100)
xlabel('\Delta O^{wet}_{2,a} / %')
ylabel('count');

subplot(122)
histogram(d_O2air_a*100)
xlabel('\Delta O^{air}_{2,a}/ %')
ylabel('count');

figure()
plot(d_O2wet_w*100+100,d_O2wet_a*100+100,'*')




