clear; close all; clc;

load('./mat_files/RGclim_labsea_012004_052020.mat');
lon = RGclim.lon;
lat = RGclim.lat;
mtime = RGclim.time;
pres = RGclim.pres;
bathy_mask = RGclim.bathy_mask;

T = RGclim.temp;  % in-situ
S = RGclim.salin; % pss-78

M = length(lat); N = length(lon); C = length(mtime); K = length(pres);

%% Apply Bathy Mask
for j = 1:C
    for i = 1:K
        temp =  T(:,:,i,j);
        salt =  S(:,:,i,j);
        temp(bathy_mask(:,:,i)~=1)=NaN;
        salt(bathy_mask(:,:,i)~=1)=NaN;
        T(:,:,i,j) = temp;
        S(:,:,i,j) = salt;
    end
end

ref_rho = 0.05;
for j = 1:M
    for i = 1:N
        for t = 1:C
            temp = T(i,j,:,t); temp = temp(:);
            salt = S(i,j,:,t); salt = salt(:);
            rho = sw_dens0(salt,temp);
            
            idx = find(~isnan(rho));
            
            if length(idx)>10
                rho(~idx)=[];
                temp(~idx)=[];
                salt(~idx)=[];
                pres(~idx)=[];
                PMLD(i,j,t)=calcmld(pres,rho,ref_rho);
            else
                PMLD(i,j,t)=NaN;
            end
        end
    end
end

[yy,mm,~] = datevec(mtime);
id = (mm ==1 |mm==2| mm == 3);


%% LOAD Sunfish
path_name = './mat_files/';
var_name  = 'sunfish_data';
load(fullfile(path_name,[var_name,'_oxy_qc.mat']));
dat = sunfish; clear sunfish
dat.gridded.rho = sw_dens(dat.gridded.salinity,dat.gridded.temperature,dat.gridded.pressure);
dat.gridded.oxygen_adjusted_umolkg = dat.gridded.oxygen_adjusted./(dat.gridded.rho/1000);
dat.gridded.oxygen_raw_umolkg = dat.gridded.oxygen_raw./(dat.gridded.rho/1000);
dat.gridded.sigma_t = sw_dens0(dat.gridded.salinity,dat.gridded.temperature)-1000;
dat1 = dat; clear dat;

path_name = './mat_files/';
var_name  = 'pearldiver_data';
load(fullfile(path_name,[var_name,'_oxy_qc.mat']));
dat = pearldiver; clear pearldiver
dat.gridded.rho = sw_dens(dat.gridded.salinity,dat.gridded.temperature,dat.gridded.pressure);
dat.gridded.oxygen_adjusted_umolkg = dat.gridded.oxygen_adjusted./(dat.gridded.rho/1000);
dat.gridded.oxygen_raw_umolkg = dat.gridded.oxygen_raw./(dat.gridded.rho/1000);
dat.gridded.sigma_t = sw_dens0(dat.gridded.salinity,dat.gridded.temperature)-1000;
dat2 = dat; clear dat;


figure(); 

t = tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

% Pearldiver
nexttile(t,1,[1 1]); hold on;
idgl = dat2.gridded.lon<-51.3 & dat2.gridded.lat>56 & dat2.gridded.timeg<datenum(2020,05,20);
dat2.gridded.dist = pathdist(dat2.gridded.lat,dat2.gridded.lon,'km');
dt = gradient(dat2.gridded.timeg(idgl))*24;
dx = gradient(dat2.gridded.dist(idgl));
[N,~,~] = histcounts2(dx,dt,[0:0.1:5],[0:0.1:3.5]);
pcolor([0.1:0.1:5]',[0.1:0.1:3.5]',N'); shading flat
plot(nanmedian(dx),nanmedian(dt),'o','MarkerFaceColor','m','MarkerEdgeColor','k','MarkerSize',10);
text(1,1,{['median dx = ',num2str(round(nanmean(dx),2)),' km'];['median dt = ',num2str(round(nanmean(dt),2)),' hrs']});
xlabel('dx / km')
ylabel('dt / hours');
colormap(gca,cmocean('tempo',10));
caxis([0 10])
cb = colorbar;
ylabel(cb,'Profile Counts')
formatplot;
xlim([0 5]);
ylim([0.1 3.5])
title('Pearldiver')



k1_mor = [-52.6412833, 56.5699833];
k78910_mor = [-51.5516833,	52.8436833;...
              -51.3082167,	52.96135; ...
              -50.8669167,	53.1339833; ...
              -50.2399167,	53.3913];

nexttile(t,2,[2 1]); hold on
contourf(lon,lat,nanmean(PMLD(:,:,id),3)');
cb1=colorbar; cb1.TickLength = 0.045;
borders('Canada','facecolor','k'); 
ylabel(cb1,'Z_{mld} / m')
colormap(gca,cmocean('amp',8));
caxis([200 1800])
formatplot
title('Glider Missions 2020 - 2022')
[lat,lon] = meshgrid(44:1/12:60.5,-66:1/12:-44);
Z = topo_interp(lat,lon);
contour(lon,lat,Z,'LevelList',[-3000:1000:-1000],'color','m')
idgl = ~isnan(dat2.longitude);
h1=plot(dat2.longitude(idgl),dat2.latitude(idgl),'-g','linewidth',2);%'o','MarkerEdgeColor','none','MarkerFaceColor','r','MarkerSize',3)
idgl = ~isnan(dat1.longitude);
h2=plot(dat1.longitude(idgl),dat1.latitude(idgl),'-b','linewidth',2);%'o','MarkerEdgeColor','none','MarkerFaceColor','b','MarkerSize',3)
h3=plot(k78910_mor(:,1),k78910_mor(:,2),'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
h4=plot(k1_mor(:,1),k1_mor(:,2),'p','MarkerEdgeColor','k','MarkerFaceColor','y','MarkerSize',10);
legend([h1 h2 h3 h4],{'Pearldiver','Sunfish','53N Array','K1 Mooring'},'Location','SW')
ylim([44 60.5])

% Sunfish
nexttile(t,3,[1 1]); hold on
idgl = dat1.gridded.timeg>datenum(2022,01,01) & dat1.gridded.timeg<datenum(2022,05,15) & ~isnan(dat1.gridded.lat);
dat1.gridded.dist = pathdist(dat1.gridded.lat(idgl),dat1.gridded.lon(idgl),'km');
dt = [0;row2col(diff(dat1.gridded.timeg(idgl)))]*24;
dx = [0;row2col(diff(dat1.gridded.dist))];
[N,Xedges,Yedges] = histcounts2(dx,dt,[0:0.1:6],[0:0.1:6]);
pcolor([0.1:0.1:6]',[0.1:0.1:6]',N'); shading flat
plot(nanmedian(dx),nanmedian(dt),'o','MarkerFaceColor','m','MarkerEdgeColor','k','MarkerSize',10);
text(1,1,{['median dx = ',num2str(round(nanmean(dx),2)),' km'];['median dt = ',num2str(round(nanmean(dt),2)),' hrs']});
xlabel('dx / km')
ylabel('dt / hours');
colormap(gca,cmocean('tempo',10));
caxis([0 10])
cb = colorbar;
% cb.Location = 'SouthOutside';
ylabel(cb,'Profile Counts')
formatplot;
xlim([0 5]);
ylim([0.1 4])
title('Sunfish')


save_figure(gcf,'./plots/deployment_map',[7.5 4],'.png','300')