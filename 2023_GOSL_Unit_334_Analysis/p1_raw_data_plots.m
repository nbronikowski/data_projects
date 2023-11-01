close all; clc;
a = 1;
% figure()
% gdat.chla(gdat.chla<0)=NaN;
% scatter(gdat.timeDateNum,-gdat.depth,20,gdat.chla,'filled');
% cb=colorbar; ylabel(cb,'Chlorophyll\_a (\mug L^{-1})');
% ylabel('depth(m)'); formatplot;
% set(gca,'XTick',datenum(2023,08,01):2:datenum(2023,09,24))
% datetick('x','dd-mmm','keepticks')
% xlim([datenum(2023,08,03) datenum(2023,09,25)])
% colormap(cmocean('algae'))
% caxis([0 2])
% title('Glider unit\_334 Wetlabs FLBB Sensor Data')
% save_figure(gcf,'chlorophyll_a',[7.5 5],'.png','300');

% figure()
% gdat.bb700(gdat.bb700<0 | gdat.bb700>1.2e-3)=NaN;
% scatter(gdat.timeDateNum,-gdat.depth,20,gdat.bb700.*1e3,'filled');
% cb=colorbar; ylabel(cb,'Backscatter 700nm (x10^{-3})');
% ylabel('depth(m)'); formatplot;
% set(gca,'XTick',datenum(2023,08,01):2:datenum(2023,09,24))
% datetick('x','dd-mmm','keepticks')
% xlim([datenum(2023,08,03) datenum(2023,09,25)])
% colormap(cmocean('matter'))
% % caxis([0 2])
% title('Glider unit\_334 Wetlabs FLBB Sensor Data')
% save_figure(gcf,'backscatter700',[7.5 5],'.png','300');



% figure()
% scatter(gdat.timeDateNum,-gdat.depth,20,gdat.temp,'filled');
% cb=colorbar; ylabel(cb,'Temperature (^oC)');
% ylabel('depth(m)'); formatplot;
% set(gca,'XTick',datenum(2023,08,01):2:datenum(2023,09,24))
% datetick('x','dd-mmm','keepticks')
% xlim([datenum(2023,08,03) datenum(2023,09,25)])
% colormap(cmocean('thermal'))
% caxis([0 20])
% title('Glider unit\_334 Seabird CTD Sensor Data')
% save_figure(gcf,'temperature',[7.5 5],'.png','150');
% 
% figure()
% gdat.salt(gdat.salt<20)=NaN;
% scatter(gdat.timeDateNum,-gdat.depth,20,gdat.salt,'filled');
% cb=colorbar; ylabel(cb,'Salinity (PSU)');
% ylabel('depth(m)'); formatplot;
% set(gca,'XTick',datenum(2023,08,01):2:datenum(2023,09,24))
% datetick('x','dd-mmm','keepticks')
% xlim([datenum(2023,08,03) datenum(2023,09,25)])
% colormap(cmocean('haline'))
% caxis([30 34.5])
% title('Glider unit\_334 Seabird CTD Sensor Data')
% save_figure(gcf,'salinity',[7.5 5],'.png','150');




figure()
gdat.u_dac(gdat.u_dac>0.5)=NaN;
gdat.v_dac(gdat.v_dac>0.5)=NaN;
plot(gdat.lon,gdat.lat,'.'); hold on; borders('Canada','facecolor','g');
s=quiver(gdat.lon,gdat.lat,gdat.u_dac,gdat.v_dac,'filled','Color','r');
s.AutoScale='off';
set(gca,'xlim',[-60.5 -57],'ylim',[49 51.5]);
xlim = get(gca,'XLim'); ylim = get(gca,'Ylim');
[xq,yq]=meshgrid(xlim(1):1/10:xlim(2),ylim(1):1/10:ylim(2));
zq = topo_interp(yq,xq);
contour(xq,yq,zq,'color','k','levellist',...
    [-1000,-300,-200,-100,-50],'showtext','on');
quiver(-60.3,49.2,0.3,0,'filled','Color','k','AutoScale','off');
formatplot
plot([-58.0663 -57.9818 -57.9818 -58.0663 -58.0663],[49.571 49.571 ...
    49.6384 49.6384 49.571],'k','LineWidth',1)
text(-60.3,49.3,0,'30 cm s^{-1}')

ax2=axes();
ax2.Position=[0.537 0.497 0.361 0.3833];
plot(ax2,gdat.lon,gdat.lat,'.'); hold on; borders('Canada','facecolor','g');
quiver(ax2,gdat.lon,gdat.lat,gdat.u_dac,gdat.v_dac,...
    'filled','Color','r','AutoScaleFactor',100);
formatplot
set(ax2,'xlim',[-58.0663 -57.9818],'ylim',[49.571 49.6384]);
set(ax2,'FontSize',8,'TickDir','in')
set(ax2,'XTickLabel','','YTickLabel','');
% [xq,yq]=meshgrid(xlim(1):1/20:xlim(2),ylim(1):1/20:ylim(2));
% zq = topo_interp(yq,xq);
% contour(ax2,xq,yq,zq,'color','k','levellist',...
%     [-200,-100,-150,-50,-20,-10],'showtext','on');
save_figure(gcf,'gliderTrajectory',[7.5 5],'.png','300');
