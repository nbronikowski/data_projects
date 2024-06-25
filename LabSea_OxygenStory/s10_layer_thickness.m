clear; close all; clc

%% Compute layer thickness and plot as well as for oxygen



path_name = './mat_files/';
var_name  = 'sunfish_data';
load(fullfile(path_name,[var_name,'_oxy_qc.mat']));
dat = sunfish; 
dat.gridded.rho = sw_dens(dat.gridded.salinity,dat.gridded.temperature,dat.gridded.pressure);
dat.gridded.oxygen_adjusted_umolkg = dat.gridded.oxygen_adjusted./(dat.gridded.rho/1000);
dat.gridded.oxygen_raw_umolkg = dat.gridded.oxygen_raw./(dat.gridded.rho/1000);
dat.gridded.sigma_t = sw_dens0(dat.gridded.salinity,dat.gridded.temperature)-1000;

% [~,mm,dd]=datevec(profs.time);

time_sc = scale_var(dat.gridded.timeg,1/4)';
time_bins = unique(time_sc);

sigma_t_sc   = scale_var(dat.gridded.sigma_t,0.001);
% sigma_t_bins = unique(sigma_t_sc(:))';
sigma_t_bins = (27.5:0.001:27.75)';

% Preallocate
sigma_layer_thick = NaN(length(sigma_t_bins),length(time_bins));
sigma_oxy_content = sigma_layer_thick;
 
for i = 1:length(time_bins)
   idxt = find(time_sc == time_bins(i));
   temp_pdens = sigma_t_sc(:,idxt);
   temp_press = dat.gridded.pressure(:,idxt);
   temp_oxy   = dat.gridded.oxygen_adjusted_umolkg(:,idxt); % micro mol / kg
   temp_pdens = temp_pdens(:);
   temp_press = temp_press(:);
   temp_oxy   = temp_oxy(:);
   for j = 1:length(sigma_t_bins)
       idxd = temp_pdens == sigma_t_bins(j);
       var_y = temp_press(idxd);
       var_x = temp_oxy(idxd);
       var_x(isnan(var_y))=[]; % Delete nans
       var_y(isnan(var_y))=[]; % Delete nans
       vol_temp = max(var_y)-min(var_y);       
       if isempty(vol_temp)
          sigma_layer_thick(j,i)=NaN;
       else
          sigma_layer_thick(j,i) = vol_temp;
       end
       if length(var_x)>3 && length(var_x)==length(var_y)
%            oxy_cont = trapz(var_y,var_x); % mols / m^2
           sigma_oxy_content(j,i) = nanmean(var_x); % mean oxygen concentration
       else
           sigma_oxy_content(j,i) = NaN;
       end
       
   end
end

close all


ptimes = [datenum(2022,01,30),datenum(2022,02,25),datenum(2022,03,31)];

figure(); 

t=tiledlayout(2,1,'TileSpacing','compact','Padding','compact')

nexttile; hold on
imagescn(time_bins,sigma_t_bins,sigma_layer_thick)
ylim([27.5   27.75]); caxis([20 200])
plot(ptimes,ptimes*0+27.75,'V','MarkerFaceColor',rgb('orange'),'MarkerEdgeColor','k');
vlines(ptimes,':k');
cb1=colorbar;
ylabel('\sigma_t / kg m^{-3} ')
formatplot; ytickformat('%.2f')
xlim([datenum(2022,01,01) datenum(2022,05,30)])
datetick('x','keeplimits')
ylabel(cb1,'Layer Thickness / m')

nexttile
imagescn(time_bins,sigma_t_bins,sigma_oxy_content)
ylim([27.5   27.75]); 
caxis([280 330])
load('odv_cmap.mat');
vlines(ptimes,':k');
colormap(gca,cmap); shading flat; cb2=colorbar;
ylabel('\sigma_t / kg m^{-3} ')
formatplot; ytickformat('%.2f')
xlim([datenum(2022,01,01) datenum(2022,05,30)])
datetick('x','keeplimits')
ylabel(cb2,'O_2 / \mumol kg^{-1}')

subtitle(t,'Sunfish \sigma_t Layer Thickness and Mean Oxygen Concentrations')

save_figure(gcf,'./plots/sunfish_sigma_t_thickness',[7.5 4.5],'.png','300');



%% T-S Profiles - pick 3 spots (1 before, 1 during, 1 after) and show Oxygen as Color
% dat.rho = sw_dens(dat.salinity_cor,dat.temperature,dat.pressure);
% dat.oxygen_adjusted_umolkg = dat.adjusted_oxygen_concentration./(dat.rho/1000);
% dat.sigma_t = sw_dens0(dat.salinity_cor,dat.temperature)-1000;

id_p1 = find(abs(dat.gridded.timeg-ptimes(1))<1/4); % 1 d
id_p2 = find(abs(dat.gridded.timeg-ptimes(2))<1/4); % 1 d
id_p3 = find(abs(dat.gridded.timeg-ptimes(3))<1/4); % 1 d

p1_S = nanmean(dat.gridded.salinity(25:1000,id_p1),2);
p1_T = nanmean(dat.gridded.temperature(25:1000,id_p1),2);
p1_O2= nanmean(dat.gridded.oxygen_adjusted_umolkg(25:1000,id_p1),2);

p2_S = nanmean(dat.gridded.salinity(25:1000,id_p2),2);
p2_T = nanmean(dat.gridded.temperature(25:1000,id_p2),2);
p2_O2= nanmean(dat.gridded.oxygen_adjusted_umolkg(25:1000,id_p2),2);

p3_S = nanmean(dat.gridded.salinity(25:1000,id_p3),2);
p3_T = nanmean(dat.gridded.temperature(25:1000,id_p3),2);
p3_O2= nanmean(dat.gridded.oxygen_adjusted_umolkg(25:1000,id_p3),2);

Pgrid = dat.gridded.pressure_grid(25:1000)';

xT = [0:0.01:7]; xS = [33:0.01:35];
[Sq,Tq] = meshgrid(xS, xT);
sigma_t_grid = sw_dens0(Sq, Tq)-1000;




figure(); 
subplot(1,2,1); hold on
plot(p1_T,-Pgrid,'color',rgb('orange'));
plot(p2_T,-Pgrid,'color',rgb('orange'),'linewidth',2);
plot(p3_T,-Pgrid,'color',rgb('orange'),'linewidth',3);
ylabel('Depth / m');
ax1 = gca; ax1.XColor=rgb('orange');
xlabel(ax1,'Temperature / C')
set(ax1, 'XAxisLocation', 'bottom');
set(ax1,'Position',[0.1018 0.1119 0.2375 0.7405])
xlim([2.7 4.55])

ax2 = axes('Position',ax1.Position, 'XAxisLocation','top', 'Color','none', ...
    'YAxisLocation','right', 'YColor','none','XColor',rgb('royal blue'));
linkaxes([ax1, ax2], 'y'); % Link the y-axis of both axes
xlabel(ax2, 'O_2 / \mumol kg^{-1}'); hold on
plot(p1_O2,-Pgrid,'color',rgb('royal blue'), 'Parent', ax2,'linewidth',1);
plot(p2_O2,-Pgrid,'color',rgb('royal blue'), 'Parent', ax2,'linewidth',2);
plot(p3_O2,-Pgrid,'color',rgb('royal blue'), 'Parent', ax2,'linewidth',3);
xlim([265 320])

subplot(1,2,2); hold on
contour(Sq,Tq,sigma_t_grid,'color','k','levelstep',0.02,'showtext','on')
id_s = 1:length(p1_S);
p1=scatter(p1_S(id_s),p1_T(id_s),10,p1_O2(id_s),'filled','Marker','o','MarkerEdgeColor','none')
p2=scatter(p2_S(id_s),p2_T(id_s),20,p2_O2(id_s),'filled','Marker','o','MarkerEdgeColor','none')
p3=scatter(p3_S(id_s),p3_T(id_s),40,p3_O2(id_s),'filled','Marker','o','MarkerEdgeColor','none')
load('odv_cmap.mat');
colormap(gca,cmap); cb=colorbar; 
ylabel(cb,'O_2 / \mumol kg^{-1}')
caxis([270,310]); 
ylabel('Temperature / C')
xlabel('Salinity / PSU')
ylim([3.4212    3.9749])
xlim([34.7980   34.8945])
set(gca,'Position',[0.4167 0.1100 0.3904 0.7405]);
save_figure(gcf,'./plots/sunfish_TS_oxygen',[7.5 4.5],'.png','300')
