clear; clc; close all

ncname = 'ERA5_uv_cloudcover.nc';
lon=double(ncread(ncname,'longitude'));
lat=double(ncread(ncname,'latitude'));
time=double(ncread(ncname,'time'))/24+datenum(1900,0,0);
tcc =double(ncread(ncname,'tcc'))*100;
u10 = double(ncread(ncname,'u10'));
v10 = double(ncread(ncname,'v10'));

[yyyy,mm,dd,HH,MM,SS]=datevec(time);
udd = unique(dd); % unique days



[xq,yq]=meshgrid(lon,lat);
[lon_z,lat_z,topo_z]=load_bathy([min(lat) max(lat)],[min(lon) max(lon)]);
topo_z_i = interp2(lon_z,lat_z,topo_z',xq,yq);
topo_z_i(topo_z_i>0)=NaN;
idnan = isnan(topo_z_i);



PanFac = 4; % we want 4 columns

% Create first figure and tiled layout for wind vectors
fig1=figure('Name', 'Wind Vectors');
t1 = tiledlayout(PanFac, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Create second figure and tiled layout for Ekman pumping velocity
fig2=figure('Name', 'Ekman Pumping Velocity');
t2 = tiledlayout(PanFac, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 20:31
    idx = dd == udd(i);
    
    % Wind vectors plot
    nexttile(t1); hold on;
    temp_u = nanmean(u10(:,:,idx), 3)';
    temp_v = nanmean(v10(:,:,idx), 3)';
    temp_u(idnan) = NaN;
    temp_v(idnan) = NaN;
    quiver(lon, lat, temp_u, temp_v, 2, 'color', 'b');
    borders('Canada', 'k', 'Linewidth', 1);
    title(datestr(datenum(2023, 08, udd(i)), 'yyyy-mmm-dd'));
    ylim([47, 55]);
    xlim([-61, -57]);
    formatplot;
    
    % Ekman pumping velocity plot
    nexttile(t2); hold on;
    [UE, VE, wE, dE] = ekman(xq, yq, temp_u, temp_v,'rho',temp_u*0+1022);
    imagescn(lon, lat, wE*1000); 
    borders('Canada', 'k', 'Linewidth', 1);
    colormap(gca,cmocean('balance'))
    title(datestr(datenum(2023, 08, udd(i)), 'yyyy-mmm-dd'));
    ylim([47, 55]);
    xlim([-61, -57]);
    formatplot;
    caxis([-1.5e-2 1.5e-2])

    if i == 22
        cb = colorbar;
        cb.Layout.Tile = 'east'; % To position the colorbar on the east side of the tiled layout
        ylabel(cb,'U_{W,Ek} (mm s^{-1})')

    end
end
title(t1,'ERA5 Daily Averaged 10m Winds');
figure(fig1)
save_figure(gcf,['./plots/Winds'],[9 12],'.png',300')

title(t2,'Vertical Ekman Transport');
figure(fig2)
save_figure(gcf,['./plots/Ekman'],[9 12],'.png',300')


for i = 20:31
    idx = dd == udd(i);
    
    ax = nexttile; hold on
    contourf(lon,lat,nanmean(tcc(:,:,idx),3)','LevelList',[0:10:100]); hold on
    borders('Canada','y','Linewidth',1)
    colormap(gca,cmocean('-rain',9))
    title(datestr(datenum(2023,08,udd(i)),'yyyy-mmm-dd'))

    ylim([47      55])
    xlim([-61 -57])
    formatplot

    if i == 22
        cb = colorbar;
        cb.Layout.Tile = 'east'; % To position the colorbar on the east side of the tiled layout
        ylabel(cb,'Total Cloud Cover (%)')

    end
    title(t,'ERA5 Total Daily Averaged Cloud Cover');
end
save_figure(gcf,['./plots/Cloud_Cover'],[9 12],'.png',300')
