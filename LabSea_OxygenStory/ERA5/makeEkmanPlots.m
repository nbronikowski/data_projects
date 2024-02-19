clear; close all; clc; clear all;
flist = dir('*.nc');
flist(2)=[];
ttime = double(ncread(flist.name,'time'))/24+datenum(1900,1,1);
lon  = double(ncread(flist.name,'longitude'));
lat  = double(ncread(flist.name,'latitude'));
v10  = double(ncread(flist.name,'v10')); % mean sea level pressure Pa
u10  = double(ncread(flist.name,'u10'));  % surface pressure Pa
sst  = double(ncread(flist.name,'sst'))+273.15;

% t1 = datenum(2022,01,20); t2 = datenum(2022,05,10); 
% qtime = t1:7:t2;
t1 = datenum(2022,02,20); t2 = datenum(2022,03,8); 
sc_fac = 1;
qtime = t1:sc_fac:t2;
tsc   = scale_var(ttime,sc_fac);
[long,latg]=meshgrid(lon,lat);
topoMask = topo_interp(latg,long);
landMask = topoMask>0;

% % extract sunfish data
path_name = './../mat_files/'
var_name  = 'sunfish_data';
load(fullfile(path_name,[var_name,'_oxy_qc.mat']));
dat = sunfish; clear sunfish


figure()
t=tiledlayout(4,4,'TileSpacing','compact','Padding','compact');

for i = 1:length(qtime)
    idt = tsc == qtime(i);
    temp_u10 = nanmean(u10(:,:,idt),3)';
    temp_v10 = nanmean(v10(:,:,idt),3)';
    [~,~,temp_wd,temp_EKd] = ekman(latg,long,temp_u10,temp_v10,'rho',1027);
    temp_EKd(landMask)=NaN;
    temp_wd(landMask)=NaN;

    glidt = dat.dateNum>qtime(i) & dat.dateNum<qtime(i)+sc_fac;

    
    nexttile; hold on
    imagescn(long,latg,temp_wd*1e3); 
    contour(long,latg,topoMask,'LevelList',[-3500:500:-1000],...
        'LineColor','k','LineStyle','-')
    quiver(long,latg,temp_u10,temp_v10,2,'Color',[.6 .6 .6],...
        'LineWidth',0.5,'AutoScale','off','MaxHeadSize',10)
    borders('Canada','k'); colormap(gca,cmocean('balance',50))
    plot(dat.longitude(glidt),dat.latitude(glidt),'.m','MarkerSize',5);
    xlim([-63 -47.25]); ylim([48.8750 60]);   
    caxis([-0.025 0.025]); 
    title([datestr(qtime(i),'dd.mmm'),' - ',datestr(qtime(i)+sc_fac,'dd.mmm')])
   
end

save_figure(gcf,'./../plots/sunfish_era5_1d_ekman',[7.5 8],'.png','300')
