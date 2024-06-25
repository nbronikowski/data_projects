close all

t3=datenum(2020,02,08);% eddy 2
t4=datenum(2020,02,10);

latlim = [56 56.8];
lonlim = [-54 -52.33];
[long,latg]=meshgrid([lonlim(1):1/12:lonlim(end)],[latlim(1):1/12:latlim(end)]);
zbathy = topo_interp(long,latg);

idEddy = find(Xq(1,:)>t3 & Xq(1,:)<t4);
idBeforeEddy = find(Xq(1,:)>t3-4 & Xq(1,:)<t4+4);


% NUM = Xq(1,idBeforeEddy);
% LAT = lat_xg(idBeforeEddy);
% LON = lon_xg(idBeforeEddy);
% U = xg_glider_U(idBeforeEddy);
% V = xg_glider_U(idBeforeEddy);
% latO = ;
% latO = ;
% eOut=eddyfit2d(NUM(:),LAT(:),LON(:),U(:),V(:),lato,lono,D,index,z0,za,zb)
% % eOUT= eddyridges(NUM(:),LAT(:),LON(:),FMAX,FMIN,P,M,RHO)

%% 
% id2 = find(l1.gl.time>t3-4 & l1.gl.time<t4+4);
% 
% 
% figure
% 
% h1=scatter3(l1.gl.lon(id2),l1.gl.lat(id2),-l1.gl.depth(id2),40,l1.gl.O2_sat(id2), ...
%    'filled');
% colormap(cmocean('oxy'))
% caxis([85 100])


%%
figure()

m_proj('mercator','latitude',latlim,'longitude',lonlim); hold on
colormap(cmocean('balance')); 
m_contour(long,latg,zbathy,'levellist',[-6000:1000:-1000,-1000:100:-10],...
       'showtext','on','color','k'); hold on
m_grid('box','on','tickdir','out','yaxislocation','left',...
       'linestyle','-','fontsize',12);

m_scatter(lon_xg(idBeforeEddy),lat_xg(idBeforeEddy),'CData',eta_steric_surf(idBeforeEddy),...
        'MarkerFaceColor','flat','SizeData',50,'Marker','s');

m_vec(0.8,lon_xg(idBeforeEddy),lat_xg(idBeforeEddy),xg_glider_U(idBeforeEddy),xg_glider_V(idBeforeEddy),'k','edgeclip','on','linewidth',0.5)

m_vec(0.8,-52.9, 56.15,0.2,0,'k')
m_text(-52.9, 56.13,'20 m s^{-1}')
title('Pearldiver Glider Track and Depth Averaged Currents')

% caxis([0.4 0.46])
cb=colorbar;
ylabel(cb,'\eta (Steric) / m')

save_figure(gcf,'map_glider_UV',[7.5 5],'.png','300')
