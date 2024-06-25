close all; clear; clc

path_name = './../mat_files/'; var_name  = 'pearldiver_data';
load(fullfile(path_name,[var_name,'_oxy_qc.mat']));
dat = pearldiver; clear pearldiver;


t3=datenum(2020,02,07,12,0,0);% eddy 2
t4=datenum(2020,02,10);


id = find(dat.dateNum>t3 & dat.dateNum<t4);

dat.SA = gsw_SA_from_SP(dat.salinity_cor,dat.pressure,dat.longitude,dat.latitude);
dat.Ctemp = gsw_CT_from_t(dat.SA,dat.temperature,dat.pressure);
dat.rho = gsw_rho(dat.SA,dat.Ctemp,dat.pressure);
dat.ptemp0 = gsw_pt0_from_t(dat.SA,dat.temperature,dat.pressure);
dat.O2mmolkg = dat.adjusted_oxygen_concentration./(dat.rho/1000);

% scatter(dat.range(id),-dat.P(id),20,dat.O2_AOU(id),'filled');

dat.range = pathdist(dat.latitude,dat.longitude,'km');
x = dat.range(id); 
x = x-x(1);
y = dat.pressure(id);
x2 = scale_var(x,1);

pg = 0:1:1020;
[Xq,Yq] = meshgrid([nanmin(x):0.25:59.3],pg);

temp_intps = barnes(x,y,dat.temperature(id),Xq,Yq,5,20,5);
salt_intps = barnes(x,y,dat.salinity_cor(id),Xq,Yq,5,20,5);
oxy_intps = barnes(x,y,dat.O2mmolkg(id),Xq,Yq,5,20,5);


pg = 0:1:1020; pg = pg(:);
[~,ux]  = pgrid_columns(x2,y,dat.temperature(id),pg);
% salt_pg  = pgrid_columns(x2,y,dat.salinity_cor(id),pg);
lon_pg = pgrid_columns(x2,y,dat.longitude(id),pg);
lat_pg = pgrid_columns(x2,y,dat.latitude(id),pg);
% [Xq,Yq] = meshgrid([nanmin(x):0.25:59.3],pg);

% interpolate to x grid + smooth 
Nr = 1; Nc = 1;
% temp_intps=intp_x_dim(ux,temp_pg,Xq,Yq,Nr,Nc);
% salt_intps=intp_x_dim(ux,salt_pg,Xq,Yq,Nr,Nc);
long_intps=intp_x_dim(ux,lon_pg,Xq,Yq,Nr,Nc);
latg_intps=intp_x_dim(ux,lat_pg,Xq,Yq,Nr,Nc);

% temp_intps = filt2(temp_intps,0.25,4.5,'lp'); % 9.5 before
% salt_intps = filt2(salt_intps,0.25,4.5,'lp');

SA_intps = gsw_SA_from_SP(salt_intps,Yq,nanmean(long_intps,1),nanmean(latg_intps,1));
CT_intps = gsw_CT_from_t(SA_intps,temp_intps,Yq);

steric_height = gsw_steric_height(SA_intps,CT_intps,Yq,1000);
eta_surf = nanmean(steric_height(1:20,:),1);

dens0 = gsw_sigma0(SA_intps,CT_intps);

MLD = mixed_layer(salt_intps,temp_intps,pg,[],0.01);





