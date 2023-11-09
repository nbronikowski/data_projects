
clear; clc;
load('pearldiver_RToxy.mat')

%% Gridding and interpolation
pg = 0:1:1025; 
[~,~,pg_time] = pgrid_columns(data.profile_index,data.pressure,data.timeDateNum,pg);
[~,~,pg_o2conc] = pgrid_columns(data.profile_index,data.pressure,data.oxygen_concentration,pg);
[~,~,pg_o2conc_RT] = pgrid_columns(data.profile_index,data.pressure,data.oxygen_concentration_corr,pg);



%     pg_o2conc_filled = inpaint_nans(pg_o2conc,3);  % Interpolate using piecewise cubic Hermite interpolating polynomial
%     pg_o2conc_filt = medfilt2(pg_o2conc_filled,[1 2],'symmetric');
%     pg_o2conc_filt(isnan(pg_o2conc)) = NaN;
% 
% %     figure(); hold on
% %     h1=plot(optode_conc_interp,-pressure_interp,'.','Color',[.6 .6 .6]);
% %     h2=plot(pg_o2conc_filt,-pg,'-k');
% %     legend({'raw sci\_ox4 (salinity compensated)','response-time corrected'},...
% %         'location','SW');
% %     formatplot
% %     ylim([-max(pg) 0])
% 
% %     ylabel('Pressure (dbar)')
% %     xlabel('Oxygen Concentration (\mumol L^{-1})')
% %     title('Aanderaa Optode 4831 w. Fast Foil Correction')
% %     save_figure(gcf,'oxy_correction',[6 5],'.png','300')
% % 
% %     pg_o2conc_filt_ds=spike_test(pg_o2conc_filt,100,3,0);
% %     figure()
% %     subplot(121)
% %     imagescn(pg_o2conc)
% %     
% %     subplot(122)
% %     imagescn(pg_o2conc_filt_ds)
% 
%     
%     gridded_oxy  = pg_o2conc_filt(:);
%     time_oxy = pg_time(:);
%     [time_oxy,idor] = sort(time_oxy);
%     gridded_oxy = gridded_oxy(idor);
%     [time_oxy,idu]=unique(time_oxy);
%     gridded_oxy = gridded_oxy(idu);
%     data.oxygen_concentration_gridded_corr(seg_idx) =interp1gap(time_oxy,gridded_oxy,timeVecSec,30);
% %     
% %     


% data.oxygen_solubility = O2solubility(data.temperature,data.salinity_corrected);
% data.oxygen_saturation = data.oxygen_concentration_gridded_corr./data.oxygen_solubility*100;
% save('glider_data_oxy_processed.mat','data');





