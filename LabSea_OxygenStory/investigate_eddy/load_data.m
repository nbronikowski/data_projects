clear ; clc; close all; clear all;

load('./../../Heat_Salt_Paper/l2_corrected_profiles.mat')
load('./../pearldiver_corrected.mat')
figure()
imagescn(gl.time,-gl.zgrid,gl.O2_mmolL);
colormap(cmocean('thermal'))
caxis([280 310])

tlims = [datenum(2020,01,15),datenum(2020,02,28)];
xlim(tlims)
datetick('x','dd-mm','keeplimits')

idt = find(gl.time>datenum(2020,01,15) & gl.time<datenum(2020,05,20));
plims = [gl.profile_index(idt(1)), gl.profile_index(idt(end))];
tlims = [gl.time(idt(1)),gl.time(plims(2))];
gl.S = sgolayfilt(gl.S,1,7,[],2);
gl.T = sgolayfilt(gl.T,1,7,[],2);
gl.depth = gl.zgrid;

gl.salt_SA = gsw_SA_from_SP(gl.S,gl.depth,gl.lon,gl.lat);
gl.ptemp0  = gsw_pt0_from_t(gl.salt_SA,gl.T,gl.depth);
gl.ctemp   = gsw_CT_from_pt(gl.salt_SA,gl.ptemp0);
gl.sigma0  = gsw_sigma0(gl.salt_SA,gl.ctemp);
gl.sigma1  = gsw_sigma1(gl.salt_SA,gl.ctemp);
gl.rho     = gsw_rho(gl.salt_SA,gl.ctemp,gl.depth);
gl.mld = mixed_layer(gl.S,gl.T,gl.zgrid,gl.sigma0,0.01);
gl.mld_smoothed = filt1('lp',...
    interp1(gl.time(~isnan(gl.mld)),gl.mld(~isnan(gl.mld)),gl.time,'linear','extrap'),...
    'fc',1/11,'fs',1/nanmedian(diff(gl.time)));
gl.eta_steric=gsw_steric_height(gl.salt_SA,gl.ctemp,gl.P,1000);

gl.O2_sol = gsw_O2sol_SP_pt(gl.S,gl.ptemp0).*gl.rho/1000;
gl.O2_AOU = gl.O2_mmolL - gl.O2_sol;

%% interpolate with x dimension (time)
% load('./../pearldiver_corrected.mat');
ux = unique(gl.time);
length(ux)
length(gl.time)

t1 = datenum(2020,01,01);
t2 = datenum(2020,06,20);
[Xq,Yq]=meshgrid(t1:1/8:t2,gl.zgrid);

Nr = 1; Nc = 0;
ctemp_xg = intp_x_dim(ux,gl.ctemp,Xq,Yq,Nr,Nc);
eta_steric_xg = intp_x_dim(ux,gl.eta_steric,Xq,Yq,Nr,Nc);
saltSa_xg = intp_x_dim(ux,gl.salt_SA,Xq,Yq,Nr,Nc);
salt_xg = intp_x_dim(ux,gl.S,Xq,Yq,Nr,Nc);
temp_xg = intp_x_dim(ux,gl.T,Xq,Yq,Nr,Nc);
theta_xg = intp_x_dim(ux,gl.ptemp0,Xq,Yq,Nr,Nc);
oxy_xg = intp_x_dim(ux,gl.O2_mmolL,Xq,Yq,Nr,Nc);
dens_xg = intp_x_dim(ux,gl.DENS0,Xq,Yq,Nr,Nc);

aou_xg = intp_x_dim(ux,gl.O2_AOU,Xq,Yq,Nr,Nc);
mld_xg = intp_x_dim(ux,repmat(gl.mld',length(gl.zgrid),1),Xq,Yq,Nr,Nc);
lon_xg = intp_x_dim(ux,repmat(gl.lon',length(gl.zgrid),1),Xq,Yq,Nr,Nc);
lat_xg = intp_x_dim(ux,repmat(gl.lat',length(gl.zgrid),1),Xq,Yq,Nr,Nc);
time_xg = intp_x_dim(ux,repmat(gl.time',length(gl.zgrid),1),Xq,Yq,Nr,Nc);


nnan_time = unique(gl.time(all(isnan(gl.O2_AOU))));
for i = 1:length(nnan_time)
    idx_gaps = find(abs(Xq(1,:)-nnan_time(i))<1/40);
    temp_xg(:,idx_gaps)=NaN;
    ctemp_xg(:,idx_gaps)=NaN;
    theta_xg(:,idx_gaps)=NaN;
    saltSa_xg(:,idx_gaps)=NaN;
    salt_xg(:,idx_gaps)=NaN;
    dens_xg(:,idx_gaps)=NaN;
    eta_steric_xg(:,idx_gaps)=NaN;
    oxy_xg(:,idx_gaps)=NaN;
    aou_xg(:,idx_gaps)=NaN;
    mld_xg(:,idx_gaps)=NaN;
end
mld_xg = nanmean(mld_xg,1);
lat_xg = nanmean(lat_xg,1);
lon_xg = nanmean(lon_xg,1);
time_xg= nanmean(time_xg,1);

eta_steric_surf = nanmean(eta_steric_xg(1:20,:),1);


%% Eddy Data
l1=load('./../../Heat_Salt_Paper/l1_corrected_vector_data.mat');
id = ~isnan(l1.gl.water_u);
gl1.u = l1.gl.water_u(id);
gl1.v = l1.gl.water_v(id);
gl1.lon = l1.gl.lon(id);
gl1.lat = l1.gl.lat(id);
gl1.time = l1.gl.time(id);


xg_glider_U = mean_interp(gl1.time,gl1.u,Xq(1,:));
xg_glider_V = mean_interp(gl1.time,gl1.v,Xq(1,:));

idnan = isnan(temp_xg(100,:));
xg_glider_U(idnan)=NaN;
xg_glider_V(idnan)=NaN;
