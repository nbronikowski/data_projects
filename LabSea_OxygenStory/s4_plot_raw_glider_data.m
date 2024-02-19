clear; clc; close all;

path_name = './mat_files/'
var_name  = 'sunfish_data'
load(fullfile(path_name,[var_name,'_oxy_qc.mat']))

sunfish.gridded.rho = sw_dens(sunfish.gridded.salinity,sunfish.gridded.temperature,sunfish.gridded.pressure);
sunfish.gridded.oxygen_adjusted_umolkg = sunfish.gridded.oxygen_adjusted./(sunfish.gridded.rho/1000);
sunfish.gridded.sigma_t = sw_dens0(sunfish.gridded.salinity,sunfish.gridded.temperature);
dat = sunfish; clear sunfish;


%% 
t1 = datenum(2022,01,01); t2 = datenum(2022,05,30); 
pg = 0:1:1030; pg = pg(:);
[Xq,Yq]=meshgrid(t1:4/24:t2,pg);

[~,idx]=deleteAlmostEmptyColumns(dat.gridded.oxygen_adjusted_umolkg,dat.gridded.pressure_grid);
OXY  = interp2(dat.gridded.timeg(idx),dat.gridded.pressure_grid,dat.gridded.oxygen_adjusted_umolkg(:,idx),Xq,Yq);
SALT = interp2(dat.gridded.timeg,dat.gridded.pressure_grid,dat.gridded.salinity,Xq,Yq);
TEMP = interp2(dat.gridded.timeg,dat.gridded.pressure_grid,dat.gridded.temperature,Xq,Yq);
SIGMA = interp2(dat.gridded.timeg,dat.gridded.pressure_grid,dat.gridded.sigma_t,Xq,Yq);

% Reintroduce gaps
timeDiff = diff(dat.gridded.timeg);
gapStartIndices = find(timeDiff > 1); % 1/4 day threshold
gapEndIndices = gapStartIndices + 1;
gapStartTimes = dat.gridded.timeg(gapStartIndices);
gapEndTimes = dat.gridded.timeg(gapEndIndices);
for i = 1:length(gapStartTimes)
    gapStart = gapStartTimes(i);
    gapEnd = gapEndTimes(i);
    gapMask = Xq > gapStart & Xq < gapEnd;
    SALT(gapMask) = NaN;
    TEMP(gapMask) = NaN;
    SIGMA(gapMask) = NaN;
end

tq = dat.gridded.timeg(idx);
timeDiff = diff(tq);
gapStartIndices = find(timeDiff > 1); % 1/4 day threshold
gapEndIndices = gapStartIndices + 1;
gapStartTimes = tq(gapStartIndices);
gapEndTimes = tq(gapEndIndices);
for i = 1:length(gapStartTimes)
    gapStart = gapStartTimes(i);
    gapEnd = gapEndTimes(i);
    gapMask = Xq > gapStart & Xq < gapEnd;
    OXY(gapMask) = NaN;
end

%%
figure()
tt=tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

nexttile;
imagescn(Xq,-Yq,TEMP);
cb1=colorbar; colormap(gca,cmocean('thermal'));
datetick('x','dd.mmm','keeplimits');
ylabel('depth')
ylabel(cb1,'T / C')
formatplot
caxis([2 4])

nexttile;
imagescn(Xq,-Yq,SALT);
cb2=colorbar; colormap(gca,cmocean('haline'));
datetick('x','dd.mmm','keeplimits');
ylabel('depth')
ylabel(cb2,'S ')
formatplot
caxis([34 34.9])

nexttile;
imagescn(Xq,-Yq,SIGMA-1000);
cb3=colorbar; colormap(gca,cmocean('dens'));
datetick('x','dd.mmm','keeplimits');
ylabel('depth')
ylabel(cb3,'\sigma_{0,t} / kg m^{-3}')
formatplot
caxis([27.5 27.85])


nexttile;
imagescn(Xq,-Yq,OXY);
cb4=colorbar;  load odv_cmap.mat; 
colormap(gca,cmap); ylabel('depth')
datetick('x','dd.mmm','keeplimits');
ylabel(cb4,'O_2 / \mumol kg^{-1}')
formatplot
caxis([280 330])

subtitle(tt,['Glider ',var_name(1:end-5),' Labrador Sea 2022 Deployment']);


save_figure(gcf,['./plots/',var_name(1:end-5),'_raw_data_gridded_timeseries'], ...
    [7.5 7],'.png','300');
