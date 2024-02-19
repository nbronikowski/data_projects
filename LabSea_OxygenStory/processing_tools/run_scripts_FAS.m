clear ; clc; close all; clear all;

% load('./../l2_corrected_profiles.mat');
% tlims = [datenum(2020,01,15),datenum(2020,05,20)];
% idt = find(gl.time>datenum(2020,01,15) & gl.time<datenum(2020,05,20));
% plims = [gl.profile_index(idt(1)), gl.profile_index(idt(end))];
% tlims = [gl.time(idt(1)),gl.time(plims(2))];
% gl.S = sgolayfilt(gl.S,1,7,[],2);
% gl.T = sgolayfilt(gl.T,1,7,[],2);
% gl.depth = gl.zgrid;
% 
% gl.salt_SA = gsw_SA_from_SP(gl.S,gl.depth,gl.lon,gl.lat);
% gl.ptemp0  = gsw_pt0_from_t(gl.salt_SA,gl.T,gl.depth);
% gl.ctemp   = gsw_CT_from_pt(gl.salt_SA,gl.ptemp0);
% gl.sigma0  = gsw_sigma0(gl.salt_SA,gl.ctemp);
% gl.sigma1  = gsw_sigma1(gl.salt_SA,gl.ctemp);
% gl.rho     = gsw_rho(gl.salt_SA,gl.ctemp,gl.depth);
% gl.mld = mixed_layer(gl.S,gl.T,gl.zgrid,gl.sigma0,0.01);
% gl.mld_smoothed = filt1('lp',...
%     interp1(gl.time(~isnan(gl.mld)),gl.mld(~isnan(gl.mld)),gl.time,'linear','extrap'),...
%     'fc',1/11,'fs',1/nanmedian(diff(gl.time)));
% gl.eta_steric=gsw_steric_height(gl.salt_SA,gl.ctemp,gl.P,1000);
% 
% gl.O2_sol = gsw_O2sol_SP_pt(gl.S,gl.ptemp0).*gl.rho/1000;
% gl.O2_AOU = gl.O2_mmolL - gl.O2_sol;

%% interpolate with x dimension (time)
load('pearldiver_corrected.mat');
ux = unique(gl.time);
t1 = datenum(2020,01,15);
t2 = datenum(2020,05,20);
[Xq,Yq]=meshgrid(t1:1/10:t2,gl.zgrid);
Nr = 1; Nc = 0;
ctemp_xg = intp_x_dim(ux,gl.ctemp,Xq,Yq,Nr,Nc);
eta_steric_xg = intp_x_dim(ux,gl.eta_steric,Xq,Yq,Nr,Nc);

salt_xg = intp_x_dim(ux,gl.S,Xq,Yq,Nr,Nc);
temp_xg = intp_x_dim(ux,gl.T,Xq,Yq,Nr,Nc);
oxy_xg = intp_x_dim(ux,gl.O2_mmolL,Xq,Yq,Nr,Nc);
aou_xg = intp_x_dim(ux,gl.O2_AOU,Xq,Yq,Nr,Nc);
mld_xg = intp_x_dim(ux,repmat(gl.mld',length(gl.zgrid),1),Xq,Yq,Nr,Nc);
lon_xg = intp_x_dim(ux,repmat(gl.lon',length(gl.zgrid),1),Xq,Yq,Nr,Nc);
lat_xg = intp_x_dim(ux,repmat(gl.lat',length(gl.zgrid),1),Xq,Yq,Nr,Nc);
time_xg = intp_x_dim(ux,repmat(gl.time',length(gl.zgrid),1),Xq,Yq,Nr,Nc);


nnan_time = unique(gl.time(all(isnan(gl.ctemp))));
for i = 1:length(nnan_time)
    idx_gaps = find(abs(Xq(1,:)-nnan_time(i))<1/12);
    temp_xg(:,idx_gaps)=NaN;
    salt_xg(:,idx_gaps)=NaN;
    eta_steric_xg(:,idx_gaps)=NaN;
    oxy_xg(:,idx_gaps)=NaN;
    aou_xg(:,idx_gaps)=NaN;
    mld_xg(:,idx_gaps)=NaN;
end
mld_xg = nanmean(mld_xg,1);
lat_xg = nanmean(lat_xg,1);
lon_xg = nanmean(lon_xg,1);
time_xg= nanmean(time_xg,1);


load('./../ERA5_pearldiver.mat');
dat.time = datenum(2020,01,15):datenum(hours(1)):datenum(2020,05,20);
dat.time = dat.time(:);
dat.temp_C = mean_interp(Xq(1,:),nanmean(temp_xg(1:10,:),1),dat.time);
dat.psal_PSU = mean_interp(Xq(1,:),nanmean(salt_xg(1:10,:),1),dat.time);
dat.dox2_umolL = mean_interp(Xq(1,:),nanmean(oxy_xg(1:10,:),1),dat.time);
dat.windspeed_ms = mean_interp(ERA5.time,sqrt(ERA5.u10.^2 + ERA5.v10.^2),dat.time);
dat.windU = mean_interp(ERA5.time,ERA5.u10,dat.time);
dat.windV = mean_interp(ERA5.time,ERA5.v10,dat.time);
dat.atmosphericpress_Pa = mean_interp(ERA5.time,ERA5.mslp,dat.time);
dat.atmosphericpress_atm = dat.atmosphericpress_Pa./101325;
dat.dox2_sol_umolkg = gsw_O2sol_SP_pt(dat.psal_PSU,dat.temp_C);


% Use David Nicholson Gas Toolbox
addpath(genpath('./gas_toolbox-master/'))

% Stanley et al. 2009
[Fd, Fc, Fp, Deq, Ks] = fas(dat.dox2_umolL*1e-3,dat.windspeed_ms,...
    dat.psal_PSU,dat.temp_C,dat.atmosphericpress_atm,'O2','S09');
Fs09 = (Fd+Fc+Fp);
Fdiff09 = Fd;
Fbubb09 = Fc+Fp;

pO2air = gasmolfract('O2').*(1-vpress(dat.psal_PSU,dat.temp_C));
pO2w   = O2ctoO2p(dat.dox2_umolL,...
    dat.temp_C,dat.psal_PSU,dat.atmosphericpress_atm)/1013.25;

% [Fd, Fc, Fp, Deq, k] = fsa(pO2w,pO2air,dat.windspeed_ms,...
%     dat.psal_PSU,dat.temp_C,'O2','L13'); %% LIANG 2013;

[F92, k] = fas_Fd(dat.dox2_umolL*1e-3,dat.windspeed_ms,...
    dat.psal_PSU,dat.temp_C,dat.atmosphericpress_atm,'O2','W92a');

% % Put back the nans
w92_O2 = mean_interp(dat.time,F92,Xq(1,:))*-86400;
s09_O2 = mean_interp(dat.time,Fs09,Xq(1,:))*-86400;
s09_O2_diff = mean_interp(dat.time,Fdiff09,Xq(1,:))*-86400;
s09_O2_bubb = mean_interp(dat.time,Fbubb09,Xq(1,:))*-86400;

windU = mean_interp(dat.time,dat.windU,Xq(1,:));
windV = mean_interp(dat.time,dat.windV,Xq(1,:));

% 
% % bar()
idnan = isnan(temp_xg(20,:));
w92_O2(idnan)=NaN;
s09_O2(idnan)=NaN;
s09_O2_diff(idnan)=NaN;
s09_O2_bubb(idnan)=NaN;

eta_steric_surf = nanmean(eta_steric_xg(1:20,:),1);

% % 
% flux_fig


% 
% 
% 
% 
% 
% 
% 
% 
% %% Net uptake calc
% % 
% % Fs09=-86400*(Fd+Fc+Fp); % mol / m2 /day
% % tq = dat.time - dat.time(1); % frac day
% % 
% % simps(tq,Fs09)
% 
% %% Compare with Glider obs directly
% 
% 
% % %% Era 5 Statistics for Uptake
% % figure()
% % subplot(121);
% % histogram(dat.windspeed_ms,0:1:26,'FaceColor','b');
% % vline(mean(dat.windspeed_ms),'--k')
% % % vline(median(dat.windspeed_ms),'--m')
% % vline(mean(dat.windspeed_ms)-1.96*std(dat.windspeed_ms),':r')
% % vline(mean(dat.windspeed_ms)+1.96*std(dat.windspeed_ms),':r')
% % title('ERA5 (JFMA) Windspeeds');
% % xlabel('u10 / m s^{-1}');
% % ylabel('Count');
% % formatplot; grid on
% % text(15,120,['mean=',num2str(mean(dat.windspeed_ms))])
% % 
% % subplot(122);
% % histogram(dat.dox2_sat-100,-15:1:10,'FaceColor','b');
% % h1=vline(mean(dat.dox2_sat-100),'--k');
% % % vline(median(dat.dox2_sat-100),'--m')
% % h2=vline(mean(dat.dox2_sat-100)-1.96*std(dat.dox2_sat-100),':r');
% % vline(mean(dat.dox2_sat-100)+1.96*std(dat.dox2_sat-100),':r')
% % title('Glider (JFMA) Supersaturation');
% % xlabel('O_2 / %');
% % % ylim([0 900])
% % formatplot; grid on
% % leg = legend([h1 h2],{'mean','95% CI'});
% % leg.Position = [0.5815 0.6996 0.1444 0.1128];
% % text(2,280,['mean=',num2str(mean(dat.dox2_sat-100))])
% % 
% % save_figure(gcf,'data_stats',[7.5 4],'.png','300')
% % 
% % 
