clear; close all; clc;

load glider_data_oxy.mat

Vlims = [-340 0];
Hlims = [nanmin(gdat.gridded.time) nanmax(gdat.gridded.time)];

figure()
VAR = gdat.gridded.temp;
ax1=subplot(311); hold on
contourf(gdat.gridded.time,-gdat.gridded.depthg,VAR,'LineStyle','none','LevelStep',0.02);
contour(gdat.gridded.time,-gdat.gridded.depthg,gdat.gridded.dens0,'ShowText','on', ...
    'LineColor','k','LevelStep',0.5)
ylim(Vlims)
xlim(Hlims)
formatplot
cb1=colorbar;
colormap(gca,cmocean('thermal'))
set(ax1,'XTickLabel','')
caxis([nanmin(VAR(:)) nanmax(VAR(:))])
ylabel('Depth (m)')
ylabel(cb1,'Temperature (^o C)')
tit=title('Glider "Batray" 2023 SWOT Calval CTD and Optode Data');

VAR = gdat.gridded.salt;
ax2=subplot(312); hold on
contourf(gdat.gridded.time,-gdat.gridded.depthg,VAR,'LineStyle','none','LevelStep',0.02);
contour(gdat.gridded.time,-gdat.gridded.depthg,gdat.gridded.dens0,'ShowText','on', ...
    'LineColor','k','LevelStep',0.5)
xlim(Hlims)
colormap(gca,cmocean('haline'))
set(ax2,'XTickLabel','')
ylim(Vlims)
formatplot
cb2=colorbar;
caxis([nanmin(VAR(:)) nanmax(VAR(:))])
ylabel('Depth (m)')
ylabel(cb2,'Salinity (PSU)')

VAR = gdat.gridded.oxy_corr;
ax3=subplot(313); hold on
load odv_cmap.mat
contourf(gdat.gridded.time,-gdat.gridded.depthg,VAR,'LineStyle','none','LevelStep',1);
contour(gdat.gridded.time,-gdat.gridded.depthg,gdat.gridded.dens0,'ShowText','on', ...
    'LineColor','k','LevelStep',0.5)
ylim(Vlims)
colormap(gca,cmap)
formatplot
cb3=colorbar;
caxis([200 nanmax(VAR(:))])
ylabel('Depth (m)')
ylabel(cb3,'Oxygen (\mumol L^{-1})')
xlim(Hlims)
datetick('x','dd.mmm','keeplimits')

ax1.Position = [0.086 0.66 0.8 0.294];
ax2.Position = [0.086 0.35 0.8 0.294];
ax3.Position = [0.086 0.04 0.8 0.294];

save_figure(gcf,['./plots/glider_CTD_O2'],[8 6],'.png','300')