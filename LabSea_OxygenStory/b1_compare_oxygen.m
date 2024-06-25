clear; clc; close all;

path_name = './mat_files/'
var_name  = 'pearldiver_data'
load(fullfile(path_name,[var_name,'_oxy_qc.mat']))
dat = pearldiver; clear pearldiver
dat.gridded.rho = sw_dens(dat.gridded.salinity,dat.gridded.temperature,dat.gridded.pressure);
dat.gridded.oxygen_adjusted_umolkg = dat.gridded.oxygen_adjusted./(dat.gridded.rho/1000);
dat.gridded.sigma_t = sw_dens0(dat.gridded.salinity,dat.gridded.temperature);

idt = dat.gridded.timeg>datenum(2020,01,01) & dat.gridded.timeg<datenum(2020,02,1);
OXY_anom = dat.gridded.oxygen_adjusted_umolkg(1:1010,:) - nanmean(dat.gridded.oxygen_adjusted_umolkg(1:1010,idt),2);


t1 = datenum(2020,01,01); t2 = datenum(2020,05,30); 
pg = 0:1:1000; pg = pg(:);
[Xq,Yq]=meshgrid(t1:2/24:t2,pg);

[OXY_anom,idx]=deleteAlmostEmptyColumns(OXY_anom,dat.gridded.pressure_grid(1:1010));
OXY_anom_intp  = interp2(dat.gridded.timeg(idx),dat.gridded.pressure_grid(1:1010),OXY_anom,Xq,Yq);

tq = dat.gridded.timeg(idx);

timeDiff = diff(tq);
gapStartIndices = find(timeDiff > 1.5); % 1/4 day threshold
gapEndIndices = gapStartIndices + 1;

% Convert indices to corresponding times
gapStartTimes = tq(gapStartIndices);
gapEndTimes = tq(gapEndIndices);

for i = 1:length(gapStartTimes)
    % Define the time range of the gap
    gapStart = gapStartTimes(i);
    gapEnd = gapEndTimes(i);

    % Find the corresponding indices in the interpolated grid
    gapMask = Xq > gapStart & Xq < gapEnd;

    % Set interpolated data to NaN within the gap
    OXY_anom_intp(gapMask) = NaN;
end



path_name = './mat_files/'
var_name  = 'sunfish_data'
load(fullfile(path_name,[var_name,'_oxy_qc.mat']))
dat1 = sunfish; clear sunfish
dat1.gridded.rho = sw_dens(dat1.gridded.salinity,dat1.gridded.temperature,dat1.gridded.pressure);
dat1.gridded.oxygen_adjusted_umolkg = dat1.gridded.oxygen_adjusted./(dat1.gridded.rho/1000);
dat1.gridded.sigma_t = sw_dens0(dat1.gridded.salinity,dat1.gridded.temperature);

t1 = datenum(2022,01,01); t2 = datenum(2022,05,30); 
pg = 0:1:1000; pg = pg(:);
[Xq1,Yq1]=meshgrid(t1:4/24:t2,pg);

idt1 = dat1.gridded.timeg>datenum(2022,01,01) & dat1.gridded.timeg<datenum(2022,02,1);
OXY_anom1 = dat1.gridded.oxygen_adjusted_umolkg(1:1010,:) - nanmean(dat1.gridded.oxygen_adjusted_umolkg(1:1010,idt1),2);

[OXY_anom1,idx]=deleteAlmostEmptyColumns(OXY_anom1,dat1.gridded.pressure_grid(1:1010));
OXY_anom_intp1  = interp2(dat1.gridded.timeg(idx),dat1.gridded.pressure_grid(1:1010),OXY_anom1,Xq1,Yq1);

tq = dat1.gridded.timeg(idx);

timeDiff = diff(tq);
gapStartIndices = find(timeDiff > 1.5); % 1/4 day threshold
gapEndIndices = gapStartIndices + 1;

% Convert indices to corresponding times
gapStartTimes = tq(gapStartIndices);
gapEndTimes = tq(gapEndIndices);

for i = 1:length(gapStartTimes)
    % Define the time range of the gap
    gapStart = gapStartTimes(i);
    gapEnd = gapEndTimes(i);

    % Find the corresponding indices in the interpolated grid
    gapMask = Xq1 > gapStart & Xq1 < gapEnd;

    % Set interpolated data to NaN within the gap
    OXY_anom_intp1(gapMask) = NaN;
end



figure()
t=tiledlayout(2,1,'TileSpacing','compact')

nexttile; hold on
imagescn(Xq,-Yq,OXY_anom_intp);
cb1=colorbar;  
title('(a) Pearldiver Oxygen Compared to Average in January  2020')
colormap(gca,cmocean('balance',20)); ylabel('Depth / m')
xlim([datenum(2020,01,15) datenum(2020,05,15)])
datetick('x','dd.mmm','keeplimits');
ylabel(cb1,'\DeltaO_2 / \mumol kg^{-1}')
formatplot
caxis([-20 20])

nexttile
imagescn(Xq1,-Yq1,OXY_anom_intp1);
cb1=colorbar;  
title('(b) Sunfish Oxygen Compared to Average in January 2022')
colormap(gca,cmocean('balance',20)); ylabel('Depth /m')
xlim([datenum(2022,01,15) datenum(2022,05,15)])
datetick('x','dd.mmm','keeplimits');
ylabel(cb1,'\DeltaO_2 / \mumol kg^{-1}')
formatplot
caxis([-20 20])

save_figure(gcf,'./plots/oxygen_comparison_both_gliders',[7.5 4],'.png','300')
