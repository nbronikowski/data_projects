clear;clc; close all; clear all;

load('./../l1_corrected_vector_data.mat')


% % Eddy 1
[~,id1] = nanmin(abs(gl.time-737810.086002885));
[~,id2] = nanmin(abs(gl.time-737815.070364463));
 

% % Eddy 2
% [~,id1] = nanmin(abs(gl.time-737828.105326902));
% [~,id2] = nanmin(abs(gl.time-737830.9851096));
 

% Eddy 3
% [~,id1] = nanmin(abs(gl.time-737920.0));
% [~,id2] = nanmin(abs(gl.time-737927.2264532));
  
sigma = gl.DENS0(id1:id2)-1000;
lon = gl.lon(id1:id2);
lat = gl.lat(id1:id2);
depth = gl.depth(id1:id2);
temp = gl.T(id1:id2);
salt = gl.S(id1:id2);
oxy  = gl.O2_sat(id1:id2);
time = gl.time(id1:id2);
u_cur = gl.water_u(id1:id2)*100; % cm/s
v_cur = gl.water_v(id1:id2)*100; % cm/s
heading = gl.heading(id1:id2);
x_comp = cosd(heading);
y_comp = sind(heading);

% temp = sw_ptmp(salt,temp,gl.P(id1:id2),0);


%% hodograph
id = ~isnan(u_cur);
U = interp1(time(id),u_cur(id),time);
id = ~isnan(v_cur);
V = interp1(time(id),v_cur(id),time);
varss = salt;
id = ~isnan(varss) & depth<20;
C = interp1(time(id),varss(id),time);

id = ~isnan(U) | ~isnan(V) | ~isnan(C);


hodograph(U(id),V(id),C(id))


xgrid = (0:1:120)';
ygrid = (0:3:1020)';
[Xq,Yq] = meshgrid(xgrid,ygrid);

% interpolate along x, y grids
[~,idlat] = nanmin(lon);
distx = lldistkm([lat,lon],[lat(idlat),lon(idlat)]);
x = scale_var(distx,median(diff(xgrid)));
y = scale_var(depth,median(diff(ygrid)));

saltp=pgrid_columns(x,depth,salt,ygrid);
tempp=pgrid_columns(x,depth,temp,ygrid);
oxyp =pgrid_columns(x,depth,oxy,ygrid);

id = saltp>nanmax(salt(:)) | saltp<nanmin(salt(:));
saltp(id)=NaN;
id = tempp>nanmax(temp(:)) | tempp<nanmin(temp(:));
tempp(id)=NaN;
id = oxyp>nanmax(oxy(:)) | oxyp<nanmin(oxy(:));
oxyp(id)=NaN;

Nr = 3; Nc = 5; 
saltx=intp_x_dim(unique(x),saltp,Xq,Yq,Nr,Nc);
tempx=intp_x_dim(unique(x),tempp,Xq,Yq,Nr,Nc);
oxyx=intp_x_dim(unique(x),oxyp,Xq,Yq,Nr,Nc);

% figure()
% pcolor(unique(x),-ygrid,tempp); shading flat; 
% 
% figure()
% pcolor(Xq,-Yq,tempx); shading flat

sc1 = 0;
sc2 = 50;
[~,idx1] = min(abs(Xq(1,:)-sc1));
[~,idx2] = min(abs(Xq(1,:)-sc2));
Xq(:,[1:idx1,idx2:end])=[];
Yq(:,[1:idx1,idx2:end])=[];
saltx(:,[1:idx1,idx2:end])=[];
tempx(:,[1:idx1,idx2:end])=[];
oxyx(:,[1:idx1,idx2:end])=[];

Rq = Xq-Xq(:,1);

Rq = [-fliplr(Rq),Rq(:,2:end)];
Yq = [fliplr(Yq),Yq(:,2:end)];
salt2 = [(saltx),fliplr(saltx(:,1:end-1))];
temp2 = [(tempx),fliplr(tempx(:,1:end-1))];
oxy2  = [(oxyx),fliplr(oxyx(:,1:end-1))];
sigma0 =sw_dens0(salt2,temp2);

id = distx<sc1 & distx>sc2;
distx2 = distx; distx2(id) = NaN;
distx2 = distx2 - sc1;

%% HATUN PAPER - BUOYANT EDDY
re = sc2; 
D = 3;
cp = 4; % kJ
r_ref = 300;
rho_ref = 1027.73;
T_ref = 3.179;
S_ref = 34.831;

for i = 1:length(Yq(:,1))
    Salt_Flux_z(i,1) = (S_ref-nanmean(salt2(i,:)))/S_ref*(re/r_ref)^2*D;
    Heat_Flux_z(i,1) = (nanmean(temp2(i,:))-T_ref)*rho_ref*cp*(re/r_ref)^2*D;
end
id = ~isnan(Heat_Flux_z);
simps(Yq(id,1),Heat_Flux_z(id))/1000
simps(Yq(id,1),Salt_Flux_z(id)*1000)


% Eddy Stats
figure()

subplot(2,2,1); hold on
contourf(Rq,-Yq,salt2,'levelstep',0.005,'linecolor','none')
[c, h]=contour(Rq,-Yq,sigma0-1000,'k','levelstep',0.05,'showtext','on');
plot(-distx2(1:100:end),-depth(1:100:end),'.k','linewidth',0.5,'MarkerSize',1);
shading flat;  
cb=colorbar; colormap(gca,cmocean('haline'))
caxis([nanmin(salt2(:))-0.001 nanmax(salt2(:))+0.001])
ylim([-1010, 0])
xlim([-sc2 0])
xlabel('Eddy Radius / km')
ylabel('Depth / m');
title('Salinity / PSU')
formatplot
set(gca,'Position',[0.1633 0.5739 0.2356 0.3412]);


subplot(2,2,2); hold on
contourf(Rq,-Yq,temp2,'levelstep',0.01,'linecolor','none')
[c, h]=contour(Rq,-Yq,sigma0-1000,'k','levelstep',0.05,'showtext','on');
shading flat;  
cb=colorbar; colormap(gca,cmocean('thermal'))
caxis([nanmin(temp(:))-0.001 nanmax(temp(:))+0.001])
ylim([-1010, 0])
xlim([-sc2 0])
xlabel('Eddy Radius / km')
title('T / deg. C')
formatplot
set(gca,'Position',[0.6166 0.5739 0.2356 0.3412]);


subplot(2,2,3)
plot(Salt_Flux_z*1000,-Yq(:,1),'b'); hold on;
plot([0 0],[-1000 0],':k');
title('Freshwater Contribution');
xlabel('[mm m^{-1}]')
ylabel('Depth / m');
formatplot; grid on
ylim([-1000, 0])

subplot(2,2,4)
plot(Heat_Flux_z,-Yq(:,1),'r'); hold on
plot([0 0],[-1000 0],':k');
title('Heat Contribution');
xlabel('[KJ m^{-3}]');
formatplot; grid on
ylim([-1000, 0])

save_figure(gcf,['pearldiver_eddy_',datestr(time(1),'mmm.dd'),'_flux'],[7.5 7],'.png','300');

% MAP
% 
% figure();
% sc = 100; 
% 
% latlim=[nanmin(lat)-0.5/12 nanmax(lat)+0.5/12];
% lonlim=[nanmin(lon)-0.5/12  nanmax(lon)+0.5/12];
% 
% 
% hold on;
% m_proj('mercator','latitude',latlim,'longitude',lonlim)
% m_etopo2('contour',[-4500:100:0],'edgecolor','k','showtext','on');
% m_line(lon,lat,'marker','.','color',[.6 .6 .6],'linewidth',1,'linestyle',':')
% 
% id = depth<10;
% m_scatter(lon(id),lat(id),60,temp(id),'filled')
% 
% colormap(cmocean('thermal',20)); cb = colorbar;
% m_vec(sc,lon,lat,u_cur,v_cur)
% m_vec(sc,-52,56.5,50,0)
% caxis([nanmin(temp(id)) nanmax(temp(id))]);
% title(['Pearldiver track ',datestr(gl.time(id1),'dd.mmm HH:MM'),...
%     ' - ',datestr(gl.time(id2),'dd.mmm HH:MM')]);
% set(gcf,'color','w');formatplot; 
% ylabel(cb,'SST (z<10m) / deg C.','FontSize',12); 
% cb.FontSize = 12;
% 
% m_grid('box','on','tickdir','out','yaxislocation','left','linestyle',':',...
%     'linewidth',0.5,'fontsize',12);
% 
% save_figure(gcf,['pearldiver_eddy_',datestr(time(1),'mmm.dd'),'_map'],[7.5 5],'.png','300');
% 
% %% Time Series
% figure()
% scatter(datetime(time,'ConvertFrom','datenum'),-depth,30,salt,'filled');
% colormap(cmocean('haline',20));
% cb= colorbar;
% formatplot;
% ylabel('Depth / m')
% ylim([-1050 0])
% ylabel(cb,'S / PSU'); grid on
% set(gca,'Position',[0.171 0.11 0.579 0.815]);
% % datetick('x','dd','keepticks')
% % xlim([time(1)-0.25 time(end)+0.25])
% cb.Position = [0.7670    0.1111    0.0617    0.8139];
% 
% save_figure(gcf,['pearldiver_eddy_',datestr(time(1),'mmm.dd'),'_profiles'],[5.5 5],'.png','300');

