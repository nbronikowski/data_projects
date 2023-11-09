clear; clc; close all;

path_name = './mat_files/'
var_name  = 'sunfish_data'
load(fullfile(path_name,[var_name,'_clean.mat']))


sunfish.gridded.sigma_t = sw_dens0(sunfish.gridded.salinity,...
                                     sunfish.gridded.temperature);
tg = nanmean(sunfish.gridded.time,1);
tg = interp1(tg(~isnan(tg)),1:length(tg))';
pg = sunfish.gridded.pressure_grid;

figure()
subplot(411);
imagescn(tg,-pg,sunfish.gridded.temperature);
cb1=colorbar; colormap(gca,cmocean('thermal'));
datetick('x','dd.mmm','keeplimits');
ylabel('depth')
ylabel(cb1,'T / C')
formatplot
title(['Raw data from ',var_name(1:end-5),' from the Labrador Sea 2022'])

subplot(412);
imagescn(tg,-pg,sunfish.gridded.salinity);
cb2=colorbar; colormap(gca,cmocean('haline'));
datetick('x','dd.mmm','keeplimits');
ylabel('depth')
ylabel(cb2,'S ')
formatplot
caxis([34 34.9])

subplot(413);
imagescn(tg,-pg,sunfish.gridded.sigma_t);
cb3=colorbar; colormap(gca,cmocean('dens'));
datetick('x','dd.mmm','keeplimits');
ylabel('depth')
ylabel(cb3,'\sigma_{0,t} / kg m^{-3}')
formatplot
caxis([1027.5 1027.85])


subplot(414);
imagescn(tg,-pg,sunfish.gridded.oxygen_raw);
cb4=colorbar;  load odv_cmap.mat; 
colormap(gca,cmap); ylabel('depth')
datetick('x','dd.mmm','keeplimits');
ylabel(cb4,'O_2 / \mumol L^{-1}')
formatplot

save_figure(gcf,['./plots/',var_name(1:end-5),'_raw_data_gridded_timeseries'], ...
    [7.5 7],'.png','300');
