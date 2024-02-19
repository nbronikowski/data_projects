clear; clc; close all;

path_name = './mat_files/';
var_name  = 'sunfish_data';
load(fullfile(path_name,[var_name,'_oxy_qc.mat']));
dat = sunfish; clear sunfish 

dat.gridded.rho = sw_dens(dat.gridded.salinity,dat.gridded.temperature,dat.gridded.pressure);
dat.gridded.oxygen_adjusted_umolkg = dat.gridded.oxygen_adjusted./(dat.gridded.rho/1000);
dat.gridded.sigma_t = sw_dens0(dat.gridded.salinity,dat.gridded.temperature)-1000;
sigma_sq = scale_var(dat.gridded.sigma_t,0.01);

nprof=length(dat.gridded.timeg);
pmld = zeros(nprof,3);
for ix = 1:nprof
    pmld(ix,:) = mld(dat.gridded.pressure_grid(:), dat.gridded.temperature(:,ix),...
        dat.gridded.salinity(:,ix), 'metric', 'threshold', ...
        'tthresh', 0.05, 'dthresh', 0.01);
end
MLD = pmld(:,3); % density threshold to estimate MLD 
MLD(MLD>999)=NaN;

% 
% figure()
% plot(dat.longitude,dat.latitude,'.')
% ylim([52 59]); xlim([-60 -49])
% hold on; 
% topoplot([-60 -49 49 55],[],-3:1:-1,'2k-')
% 

dist21000m = lldistkm([dat.latitude(:),dat.longitude(:)],[53.333,-52.0833]);

%% Gridd data to 1 hr and load era5 data
% MEDIAN RESOLUTION 4 hrs 


idt = dat.gridded.timeg>datenum(2022,01,01) & dat.gridded.timeg>datenum(2022,01,31);
oxy_anom = dat.gridded.oxygen_adjusted_umolkg-...
    nanmean(dat.gridded.oxygen_adjusted_umolkg(:,idt),2);


oxy = nan*dat.gridded.timeg;
temp = nan*dat.gridded.timeg;
salt = nan*dat.gridded.timeg;
dens = nan*dat.gridded.timeg;
for i = 1:length(dat.gridded.timeg)
%     id =  sigma_sq(:,i) == 27.68;
    id = abs(dat.gridded.pressure(:,i)-900)<5;
    oxy(i) = nanmean(oxy_anom(id,i));
    temp(i)= nanmean(dat.gridded.temperature(id,i));
    salt(i)= nanmean(dat.gridded.salinity(id,i));
%     depth(i)=nanmean(dat.gridded.pressure(id,i));
    dens(i)=nanmean(dat.gridded.sigma_t(id,i));
end
dfrac = 1/6;


t1 = datenum(2022,01,01); t2 = datenum(2022,06,01); 

map.time  = (t1:dfrac:t2)';
map.oxy = mean_interp(dat.gridded.timeg  ,oxy,map.time,0);
map.temp = mean_interp(dat.gridded.timeg  ,temp,map.time,0);
map.salt = mean_interp(dat.gridded.timeg  ,salt,map.time,0);
% map.depth=mean_interp(dat.gridded.timeg  ,depth,map.time,0);
map.lon = mean_interp(dat.gridded.timeg  ,dat.gridded.lon(:),map.time,0);
map.lat = mean_interp(dat.gridded.timeg  ,dat.gridded.lat(:),map.time,0);
map.u   = mean_interp(dat.dateNum,dat.u,map.time,0);
map.v   = mean_interp(dat.dateNum,dat.v,map.time,0);
map.dist2_1000m = mean_interp(dat.dateNum,dist21000m,map.time,0);
map.mld = mean_interp(dat.gridded.timeg,MLD,map.time,0);

load('./mat_files/sunfish_era5_data.mat');
b=load('./mat_files/sunfish_era5_wind_ekman_data.mat');
era5.u10 = b.era5.u10;
era5.v10 = b.era5.v10;
% 
map.u10    = mean_interp(era5.time,era5.u10,map.time,0);
map.v10    = mean_interp(era5.time,era5.v10,map.time,0);
map.ekd   = mean_interp(b.era5.time,b.era5.ekd,map.time,0);
map.wd    = mean_interp(b.era5.time,b.era5.wd,map.time,0);


t3 = datenum(2022,02,01); t4 = datenum(2022,05,01); 
idt = map.time>t3 & map.time<t4;

x = map.dist2_1000m(idt); y = map.oxy(idt);
t = map.time(idt);
idd = isnan(x) | isnan(y);
x(idd)=[]; y(idd)=[];
t(idd)=[];


ft = fittype('poly1'); 
fo = fitoptions('Method', 'LinearLeastSquares', 'Robust', 'bisquare');

idnan = isnan(x) | isnan(y);
[fitresult, gof] = fit(x(~idnan),y(~idnan), ft, fo);

% Get the fitted line
x_fit = linspace(min(x)-10,max(x)+10, 200); % Generating linear space for fit line
y_fit = feval(fitresult, x_fit);

% Plotting the confidence intervals
conf_int = predint(fitresult,x_fit,0.95,'functional','on');
coeffs = coeffvalues(fitresult);

% Plot original data
% figure(); hold on;
% 
% h1=plot(x,y,'ok','MarkerFaceColor','k'); 
% h2=plot(x_fit, y_fit, '-r'); % Fitted line
% h4=plot(xlim,fliplr(ylim),'-m');
% h3=plot(x_fit, conf_int, 'r:'); % 95% prediction bounds
% % Annotating the plot with fit information
% txt = sprintf('Slope: %.3f\nIntercept: %.3f\nR^2: %.3f', coeffs(1), coeffs(2), gof.rsquare);
% text(mean(x), mean(y), txt); % Adjust the position as needed
% ylabel('O_2_{(\sigma_t=27.68 kg m3^{-3})} / \mumol kg^{-1}')
% xlabel('Distance to 1000m isobar / km')
% 
% save_figure(gcf,'./plots/sunfish_oxygen_vs_distance',[5 5],'.png','300')


%% 
% - Add cbdate to scatter plot
% change scales, tidy up and save fig.
% do same for PearlDiver
% Then next steps - why 900 m? 
%


figure()

tt=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

t3 = datenum(2022,02,01); t4 = datenum(2022,05,31); 


nexttile; hold on
h1=scatter(map.time,map.dist2_1000m,20,map.oxy,'filled');
cb1=colorbar; ylabel(cb1,'\Delta O_2 / \mumol kg^{-1}')
caxis([-20 20])
colormap(gca,cmocean('balance'))
h3=hlines(nanmean(map.dist2_1000m(idt)),':k');
xlim([datenum(2022,01,20) t4]); 
% ylim([nanmin(map.dist2_1000m(idt)) nanmax(map.dist2_1000m(idt))])
ylim([70 165])
ylabel('Track Distance to z=1000m / km');
datetick('x','dd.mmm','keeplimits');
title('(a)')


nexttile; hold 
h1=scatter(x,y,20,t,'filled'); cb2=colorbar; grid on
ylim([-20 20])
xlim([50 165])
cbdate(t3:14:t4,'dd-mmm')
h2=plot(x_fit, y_fit, '-r','LineWidth',2); % Fitted line
h4=plot(xlim,fliplr(ylim),'-m');
h3=plot(x_fit, conf_int, 'r:'); % 95% prediction bounds
txt = sprintf('Slope: %.3f\nIntercept: %.3f\nR^2: %.3f', coeffs(1), coeffs(2), gof.rsquare);
text(120, 10, txt); % Adjust the position as needed
ylabel('\Delta O_2 / \mumol kg^{-1}')
xlabel('Track Distance to z=1000m / km')
title('(b)')
legend([h2 h4],{'Robust Fit','1:1 fit'},'Location','SW');
subtitle(tt,'Glider Sunfish Oxygen Anomaly to January at 900 m')

save_figure(gcf,'./plots/sunfish_oxygen_vs_distance',[7.5 4],'.png','300')

% nexttile; 
% 
% hold on;
% y = zeros(size(map.time)); % Constant baseline
% quiver(map.time, y, -map.u, map.v, 0.5,'Color','b','AutoScale','off'); % Use scaled vectors
% xlim([t3 t4]);
% datetick('x', 'dd.mmm', 'keeplimits');
% set(gca, 'YTick', '');
% ylabel('Glider DAC');
% title('(a)');

% nexttile; hold on
% h1=scatter(map.time,map.dist2_1000m,20,map.wd*1e7,'filled');
% cb=colorbar; ylabel(cb,'W_{ek} / m s^{-1} x 10^7 ')
% caxis([-100 100])
% colormap(gca,cmocean('balance'))
% h3=hlines(nanmean(map.dist2_1000m(idt)),':k');
% xlim([t3 t4]); 
% ylabel('X_{(z=1000m)} / km');
% datetick('x','dd.mmm','keeplimits');
% title('(b)')



% nexttile; hold on
% h1=scatter(map.time,map.dist2_1000m,20,map.mld,'filled');
% cb3=colorbar; ylabel(cb3,'MLD / m')
% caxis([10 200])
% h3=hlines(nanmean(map.dist2_1000m(idt)),':k');
% xlim([t3 t4]); 
% ylabel('X_{(z=1000m)} / km');
% datetick('x','dd.mmm','keeplimits');
% title('(d)')

% subtitle(t,'Sunfish Depth Averaged Currents, Ekman, Oxygen and MLD')
% save_figure(gcf,'./plots/sunfish_april_oxygen_uv',[7.5 5],'.png','300')





% h1=plot(map.time,map.oxy,'.-','color',rgb('purple'));
% h3=hlines(nanmean(map.oxy(idt)),':k');
% h2=plot(map.time,medfilt1(map.oxy,1/dfrac),'-','color',rgb('purple'),'LineWidth',2);
% xlim([t3 t4])
% ylabel('O_2 / \mumol kg^{-1}');
% datetick('x','dd.mmm','keeplimits');
% title('(c)')
% ylim([nanmin(map.oxy(idt)) nanmax(map.oxy(idt))])

% nexttile; hold on
% h1=plot(map.time,map.salt,'.-','color',rgb('red'));
% h3=hlines(nanmean(map.salt(idt)),':k');
% h2=plot(map.time,medfilt1(map.salt,1/dfrac),'-','color',rgb('red'),'LineWidth',2);
% xlim([t3 t4])
% ylabel('S / PSU');
% datetick('x','dd.mmm','keeplimits');
% title('(d)')
% ylim([nanmin(map.salt(idt)) nanmax(map.salt(idt))])

% 
% nexttile; hold on
% h1=plot(map.time,-map.depth2771,'.-','color',rgb('black'));
% h3=hlines(-nanmean(map.depth2771(idt)),':k');
% % h2=plot(map.time,medfilt1(map.salt,12),'-','color',rgb('red'),'LineWidth',2);
% xlim([t3 t4])
% ylabel('z_{(27.71 kg m^{-3})} / m');
% datetick('x','dd.mmm','keeplimits');
% title('(e)')
% % ylim([-nanmax(map.depth2771(idt)) -nanmin(map.depth2771(idt))])
% % legend([h1, h2],{'2-hr glider data';'1-d median filter'},'Location','SW','Orientation','horizontal')
% 


% %% Correlation
% rho = 1.225; % Air density (kg/m^3) at sea level
% Cd = 1.3e-3; % Drag coefficient, varies depending on the conditions
% tau_x = rho * Cd .* map.u10 .* sqrt(map.u10.^2 + map.v10.^2);
% tau_y = rho * Cd .* map.v10 .* sqrt(map.u10.^2 + map.v10.^2);
% 
% T = 4/dfrac;
% ttime = map.time-map.time(1);
% timeDays = max(scale_var(ttime,1/dfrac));
% iTime = map.time(1):1:map.time(end)-1;
% for t = 1:timeDays %-4+1
%     if t-1 <= timeDays
%             % Select the 4-day window for tau_x and tau_y
%             window_tau_x = tau_x(t:t+T-1);
%             window_tau_y = tau_y(t:t+T-1);
%             
%             % Trapezoidal integration within the 4-day window for tau_x and tau_y
%             idnan = ~isnan(window_tau_x);
%             integral_tau_x = trapz(window_tau_x(idnan)) * dfrac*86400;
%             integral_tau_y = trapz(window_tau_y(idnan)) * dfrac*86400;
%             
%             % Sum of integrals for tau_x and tau_y is the wind impulse for this window and spatial point
%             windImp(t) = (integral_tau_x + integral_tau_y);
%             meanOxy(t) = nanmean(map.oxy(t:t+T-1));
%     end
% end
% 
% 
% figure()
% t2 = tiledlayout(2,1,"TileSpacing","compact");
% 
% idt2 = iTime>t3 & iTime<t4;
% nexttile; hold on
% plot(iTime,windImp,'-')
% % xlim([t3 t4])
% datetick('x','dd.mmm','keeplimits');
% 
% 
% nexttile; hold on
% plot(iTime,meanOxy,'-')
% % xlim([t3 t4])
% datetick('x','dd.mmm','keeplimits');
% 
% 
% % en
% % 
% % oxy_c = medfilt1(map.oxy,6); %-nanmean(map.oxy(idt));
% % dekd_c = medfilt1(map.ekd,6) % ekman depth and oxygen at surface show some agreement
% % id = isnan(oxy_c) | isnan(dekd_c);
% % oxy_c(id)=[];
% % dekd_c(id)=[];
% % 
% % 
% figure()
% plot(windImp(idt2)-nanmean(windImp(idt2)),meanOxy(idt2)-nanmean(meanOxy(idt2)),'.')
