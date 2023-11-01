close all
[z_lon, z_lat,z_topo] = load_bathy([49 51], [-60.5 -57.5]);

idNot = find(gdato.timeDateNum<timeg(1) | gdato.timeDateNum>timeg(end));
u_i = mean_interp(gdato.timeDateNum,gdato.u_dac,timeg,1);
v_i = mean_interp(gdato.timeDateNum,gdato.v_dac,timeg,1);

% 
% figure(); hold on
% contour(z_lon,z_lat,z_topo','levellist',[-500,-300,-250,-200,-150,-100,-50,-25],'color',[.4 .4 .4],'showtext','on')
% 
% plot(gdato.lon(idNot),gdato.lat(idNot),'.b')
% borders('Canada','color','k')
% ax1=scatter(long,latg,70,intg_chla,'filled','Marker','square')
% quiver(long,latg,u_i',v_i','off','color','k')
% quiver(mean(long)-6.5*1/12,mean(latg),-0.1,0,'off','color','k')
% text(mean(long)-8*1/12,mean(latg)+0.4/12,'10 cm/s')
% caxis([nanmin(intg_chla) nanmax(intg_chla)])
% cb1=colormap(cmocean('algae'));
% ylim([49.5 50.2])
% xlim([-60 -58.5])
% ytickformat('%.1f')
% xtickformat('%.1f')
% cb=colorbar;
% ylabel(cb,'integrated Chlorophyll z_{mld} - z_{mld+50} (mg m^{-2})')
% formatplot
% title('Glider Vertically Integrated Chlorophyll for August 20 - 29, 2023')
% save_figure(gcf,['./plots/glider_chlor_vs_position'],[7.5 5],'.png','300')


%% Next lets look at backscatter

[xq,yq]=meshgrid(timeg,pg);
gl_mean_lon = nanmean(pg_lon,1);
gl_mean_lat = nanmean(pg_lat,1);

for i =1:length(gl_mean_lon)
    [~,idLon]=nanmin(abs(z_lon-gl_mean_lon(i)));
    [~,idLat]=nanmin(abs(z_lat-gl_mean_lat(i)));
    glZ(i)=nanmean(z_topo(idLon,idLat));
end

id_close_floor=abs(yq+glZ)<50;
z_mean_bbp=NaN*timeg;
for i = 1:length(timeg)
    idnan = ~isnan(pg_bbp700(:,i));
    
    [max_z,id_max_z]=nanmax(pg(idnan));
    [~,id_start] = nanmin(abs(pg(idnan)-(max_obs_z-20)));
    disp(['Averaging from ',num2str(pg(id_start)),' to ',num2str(pg(id_max_z))])
    z_mean_bbp(i) = nanmax(movmedian(pg_bbp700(id_start:id_max_z,i),5));


end

figure(); hold on
contour(z_lon,z_lat,z_topo','levellist',[-500,-300,-250,-200,-150,-100,-50,-25],'color',[.4 .4 .4],'showtext','on')
plot(gdato.lon(idNot),gdato.lat(idNot),'.b')
borders('Canada','color','k')
ax1=scatter(gl_mean_lon,gl_mean_lat,70,z_mean_bbp,'filled','Marker','square')
caxis([nanmin(z_mean_bbp) nanmax(z_mean_bbp)])
colormap(gca,cmocean('matter'))
ylim([49.5 50.2])
xlim([-60 -58.5])
ytickformat('%.1f')
xtickformat('%.1f')
% caxis([0 0.02])
cb2=colorbar;
ylabel(cb2,'b_{bp,700nm} (10^{-4} m^{-1})')
formatplot
title('Max backscatter within 20m of dive inflections, Aug. 20-29, 2023')
save_figure(gcf,['./plots/glider_bbp_vs_position'],[7.5 5],'.png','300')