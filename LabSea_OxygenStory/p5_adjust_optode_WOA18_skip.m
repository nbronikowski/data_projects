clear; clc; close all;

path_name = './mat_files/'
var_name  = 'pearldiver_data'
load(fullfile(path_name,[var_name,'_clean.mat']))
dat = pearldiver;  clear pearldiver

dat.gridded.time = mean(dat.gridded.time,1,'omitnan');
dat.gridded.lon = interp1(datenum(dat.time),dat.longitude,dat.gridded.time);
dat.gridded.lat = interp1(datenum(dat.time),dat.latitude,dat.gridded.time);

% Load WOA data because we don't have a lot of data right at the surface
% and this way we can use the top 20 m for a gain correction as per ARGO QC

% load woa2018 data and compute gains ...
load('./WOA2018/woa_oxy_sat_2018.mat')

%% PROFILE BY PROFILE APPROACH
[~,mm,~]=datevec(dat.gridded.time);
dat.gridded.O2sat = O2ctoO2s(dat.gridded.oxygen_raw,dat.gridded.temperature,...
                        dat.gridded.salinity,dat.gridded.pressure);

for i=1:length(dat.gridded.profile_index)
    [~,woa_lon_idx] = min(abs(woa2018.lon-dat.gridded.lon(i)));
    [~,woa_lat_idx] = min(abs(woa2018.lat-dat.gridded.lat(i)));
    g.oxy_sat_ref(i) = woa2018.obj_oxy_sat(woa_lon_idx,woa_lat_idx,mm(i))
    g.oxy_sat_mn(i) = mean(dat.gridded.O2sat(1:20,i),'omitnan');
    g.oxy_sat_std(i) = std(dat.gridded.O2sat(1:20,i),[],'omitnan');
    g.gcoef(i) = g.oxy_sat_ref(i)/g.oxy_sat_mn(i);
end
    

O2_gains= g.gcoef; 

figure(); hold on
h1=plot(dat.gridded.time,O2_gains,'.'); 
h2=plot(dat.gridded.time,dat.gridded.time*0+nanmedian(O2_gains),'-b');
h3=plot(dat.gridded.time,dat.gridded.time*0+nanmean(O2_gains),'-r');
h4=plot(dat.gridded.time,dat.gridded.time*0+nanmean(O2_gains)-nanstd(O2_gains),'--k')
plot(dat.gridded.time,dat.gridded.time*0+nanmean(O2_gains)+nanstd(O2_gains),'--k')
legend([h1 h2 h3 h4],{'gains','median','mean','\pm 1\sigma'},'Location','best')
title('Pearldiver Glider and WOA2018 Surface O_2 Saturation Comparison')
ylabel('gain factor O^{WOA}_2 / O^{glider}_2 ' )
datetick('x','dd-mmm-yy')
save_figure(gcf,['./plots/pearldiver_woa2018_optode_gains'],[7.5 5],['.png'],'300')



% figure();
% 
% subplot(211); hold on; formatplot
% plot(dat.gridded.time,g.oxy_sat_mn,'*b','DisplayName','Glider Data')
% plot(dat.gridded.time,g.oxy_sat_ref,'ok','DisplayName','WOA2018')
% ylabel('O_2,gl / %'); 
% % xlabel('O_2,WOA2018 / %')
% datetick('x','dd-mmm-yyyy','keeplimits')
% title('Pearldiver Glider and WOA 2018 Oxygen Comparison')
% 
% subplot(212)
% plot(dat.gridded.time,g.gcoef,'*b','DisplayName','gains')
% legend; ylabel('gain coefficients'); datetick('x','dd-mmm-yyyy','keeplimits')
% 
% save_figure(gcf,['./plots/pearldiver_oxy_gains_profile'],[7.5 5],'.png','300')


% APPLY CORRECTION
% dat.gridded.O2con_adj = dat.gridded.oxygen_raw.*g.gcoef;
% dat.gridded.O2sat_adj = O2ctoO2s(dat.gridded.O2con_adj,dat.gridded.temperature,...
%                         dat.gridded.salinity,dat.gridded.pressure);

% plot(dat.gridded.O2sat_adj-dat.gridded.O2sat,'k'); hold on
% plot(,'color',[.6 .6 .6])



%% MONTHLY APPROACH
% dat.O2sat = O2ctoO2s(dat.raw_oxygen_concentration,dat.temperature,...
%                         dat.salinity,dat.pressure);
% [~,mm,~]=datevec(dat.time);
% clear g
% for i = 1:5
%     tidx = mm == i & dat.depth<20;
%     g.oxy_sat_mn(i) = mean(dat.O2sat(tidx),'omitnan');
%     g.oxy_sat_std(i) = std(dat.O2sat(tidx),[],'omitnan');
% %     plot(dat.longitude(tidx),dat.latitude(tidx))
%     g.lon(i) = mean(dat.longitude(tidx),'omitnan');
%     g.lat(i) = mean(dat.latitude(tidx),'omitnan');
%     
%     [~,woa_lon_idx] = min(abs(woa2018.lon-g.lon(i)));
%     [~,woa_lat_idx] = min(abs(woa2018.lat-g.lat(i)));
%     g.oxy_sat_ref(i)= woa2018.obj_oxy_sat(woa_lon_idx,woa_lat_idx,i);
%     g.gcoef(i) = g.oxy_sat_ref(i)/g.oxy_sat_mn(i);
%     g.time(i) = datenum(2022,i,15);
% end
% 
% figure();
% 
% subplot(211); hold on; formatplot
% plot(g.oxy_sat_ref,g.oxy_sat_mn,'*b');
% % plot(g.time,g.oxy_sat_ref,'ok','DisplayName','WOA2018')
% % leg=legend; leg.Location='best';
% ylabel('O_2,meas / %'); 
% xlabel('O_2,WOA2018 / %'); 
% % datetick('x','dd-mmm-yyyy','keeplimits')
% title('Sunfish Glider and WOA 2018 Oxygen Comparison')
% 
% subplot(212)
% plot(g.time,g.gcoef,'*b','DisplayName','gains')
% legend; ylabel('gain coefficients'); datetick('x','dd-mmm-yyyy','keeplimits')
% 
% save_figure(gcf,['./plots/sunfish_oxy_gains_monthly'],[7.5 5],'.png','300')

% APPLY CORRECTION
% dat.O2con_adj = dat.raw_oxygen_concentration.*interp1(g.time,g.gcoef,datenum(dat.time));
% dat.O2sat_adj = O2ctoO2s(dat.O2con_adj,dat.temperature,...
%                         dat.salinity,dat.pressure);
% 
% plot(dat.time,dat.O2sat,'.'); hold on
% plot(dat.time,dat.O2sat_adj,'.')
