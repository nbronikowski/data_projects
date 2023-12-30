clear; clc; close all;

path_name = './mat_files/'
var_name  = 'sunfish_data'
load(fullfile(path_name,[var_name,'_oxy_qc.mat']))
dat = sunfish; 
dat.gridded.pressure_grid=ones(size(dat.gridded.pressure)).*dat.gridded.pressure_grid';

dat.gridded.rho = sw_dens(dat.gridded.salinity,dat.gridded.temperature,dat.gridded.pressure);
dat.gridded.oxygen_adjusted_umolkg = dat.gridded.oxygen_adjusted./(dat.gridded.rho/1000);
dat.gridded.oxygen_raw_umolkg = dat.gridded.oxygen_raw./(dat.gridded.rho/1000);

dat.gridded.sigma_t = sw_dens0(dat.gridded.salinity,dat.gridded.temperature);


%% Load Argo Data
load('./ARGO/argo_oxygen.mat')
argo.sigma_t = sw_dens0(argo.PSAL,argo.TEMP);

% Co-locate glider and Argo measurements and store result
idx = find(argo.TIME>=dat.gridded.timeg(1) & argo.TIME<=datenum(2022,05,20));

GL_VAR_LIST = fieldnames(dat.gridded);
AR_VAR_LIST = fieldnames(argo);
argo.pgrid = ones(size(argo.PRES)).*argo.pgrid;
m= 1;
for i = 1:length(argo.TIME(idx))
    VarCount = find(~isnan(argo.TEMP(:,idx(i))));
    if ~isempty(VarCount)
        dkm=lldistkm([dat.gridded.lat(:),dat.gridded.lon(:)],...
            [argo.LATITUDE(idx(i)),argo.LONGITUDE(idx(i))]);
        dtime = abs(dat.gridded.timeg - argo.TIME(idx(i)));
        dkm=dkm(:); dtime = dtime(:);
        idOutside = dtime>30 | dkm>30;
        dtime(idOutside)=NaN;
        dkm(idOutside)=NaN;
        [val,glIDX] = min(dkm(:)+dtime(:),[],'omitnan');
        if ~isempty(glIDX) && ~isnan(val)
            val
            dkm(glIDX)
            dtime(glIDX)
            for j = 1:length(GL_VAR_LIST)
                glProfs.(GL_VAR_LIST{j})(:,m)=dat.gridded.(GL_VAR_LIST{j})(:,glIDX);
            end
            for k = 1:length(AR_VAR_LIST)
                glProfs.argo.(AR_VAR_LIST{k})(:,m) = argo.(AR_VAR_LIST{k})(:,idx(i));
            end
            m = m+1;
        end
    end
end



%% Grid in density space
density_grid = 1027.6:0.02:1027.75;
gl_dens_sc=scale_var(glProfs.sigma_t,0.02);
ar_dens_sc=scale_var(glProfs.argo.sigma_t,0.02);

arg_oxy_dens_grd = NaN(length(density_grid),min(size(glProfs.argo.DOX2_ADJUSTED)));
gld_oxy_dens_grd = arg_oxy_dens_grd;
for i=1:length(density_grid)

    for j = 1:min(size(glProfs.argo.DOX2_ADJUSTED))
        
        id_arg = ar_dens_sc(:,j)==density_grid(i) & abs(glProfs.sigma_t(:,j)-glProfs.argo.sigma_t(:,j))<0.02;
        arg_oxy_dens_grd(i,j)=nanmean(glProfs.argo.DOX2_ADJUSTED(id_arg,j));
        arg_oxy_dens_std(i,j)=nanstd(glProfs.argo.DOX2_ADJUSTED(id_arg,j));
        id_gld = gl_dens_sc(:,j)==density_grid(i) & abs(glProfs.sigma_t(:,j)-glProfs.argo.sigma_t(:,j))<0.02;
        gld_oxy_dens_grd(i,j)=nanmean(glProfs.oxygen_adjusted_umolkg(id_gld,j));
        gld_oxy_dens_std(i,j)=nanstd(glProfs.oxygen_adjusted_umolkg(id_gld,j));
    end
end




ft = fittype('poly1'); 
fo = fitoptions('Method', 'LinearLeastSquares', 'Robust', 'bisquare');


x = arg_oxy_dens_grd(:); y = gld_oxy_dens_grd(:);
x_std = arg_oxy_dens_std(:); y_std = gld_oxy_dens_std(:);
idnan = isnan(x) | isnan(y);

[fitresult, gof] = fit(x(~idnan),y(~idnan), ft, fo);

% Get the fitted line
x_fit = linspace(250,310, 200); % Generating linear space for fit line
y_fit = feval(fitresult, x_fit);

% Plotting the confidence intervals
conf_int = predint(fitresult,x_fit,0.95,'functional','on');
coeffs = coeffvalues(fitresult);


% Plot original data
figure(); hold on;

h1=plot(x,y,'ok','MarkerFaceColor','k'); 
errorbar(x, y, x_std, 'horizontal', 'LineStyle', 'none', 'Color', 'k');
errorbar(x, y, y_std, 'vertical', 'LineStyle', 'none', 'Color', 'k');
h2=plot(x_fit, y_fit, '-r'); % Fitted line
h4=plot(xlim,ylim,'-m');
h3=plot(x_fit, conf_int, 'r:'); % 95% prediction bounds
% Annotating the plot with fit information
txt = sprintf('Slope: %.3f\nIntercept: %.3f\nR^2: %.3f', coeffs(1), coeffs(2), gof.rsquare);
text(295, 290, txt); % Adjust the position as needed

% Labels and title (customize as needed)
xlabel('O^{Argo}_2 / \mumol kg^{-1}');
ylabel('O^{Glider}_2 / \mumol kg^{-1}');
title('Sunfish vs Argo Oxygen in same \sigma_t bins (0.02 kg m^{-3})');

legend([h1; h2; h4 ;h3(1)],{'Data';'Robust Fit';'1:1 Fit';'95% FIT CI'},'Location','best')
grid on; formatplot

save_figure(gcf,['./plots/Sunfish_vs_Argo_Oxygen'],[5.5 4],'.png','300')



figure()
plot(glProfs.sigma_t-glProfs.argo.sigma_t,glProfs.pressure)



%% SAVE Data and make a residual plot for
% - glider vs argo oxygen and residuals over density space
% - and profiles vs depth




%% PLOT RESULTING STRUCTURE


tARGO = ones(size(glProfs.argo.TEMP)).*glProfs.argo.TIME;


figure()
t = tiledlayout(2,2);


nexttile; hold on
colMAP = flipud(seminfhaxby(length(glProfs.timeg)+2));
for i = 1:length(glProfs.timeg)
    plot(glProfs.salinity(:,i),-glProfs.pressure(:,i),'Color',colMAP(i,:),'LineWidth',1.5);
    plot(glProfs.argo.PSAL_ADJUSTED(:,i),-glProfs.argo.pgrid(:,i),'Color',colMAP(i,:));
end
ylim([-1000 0]); formatplot
xlabel('S')
title('(a)');

nexttile; hold on

for i = 1:length(glProfs.timeg)
    h1=plot(glProfs.temperature(:,i),-glProfs.pressure(:,i),'Color',colMAP(i,:),'LineWidth',1.5);
    h2=plot(glProfs.argo.TEMP_ADJUSTED(:,i),-glProfs.argo.pgrid(:,i),'Color',colMAP(i,:));
end
ylim([-1000 0]); formatplot
legend([h1 h2],{'Glider','Argo'},'Location','best');
xlabel('T / C')
formatplot
title('(b)')

nexttile; hold on

for i = 1:length(glProfs.timeg)
    plot(glProfs.oxygen_adjusted_umolkg(:,i),-glProfs.pressure(:,i),'Color',colMAP(i,:),'LineWidth',1.5);
    plot(glProfs.argo.DOX2_ADJUSTED(:,i),-glProfs.argo.pgrid(:,i),'Color',colMAP(i,:));
end
ylim([-1000 0]); formatplot
xlabel('O_2 / \mumol kg^{-1}')
formatplot
title('(c)')

nexttile; hold on
topoplot([-51.5 -49.5 50.7   54.9],[],[0:-1/2:-5],'-k')
plot(dat.gridded.lon,dat.gridded.lat,'.','MarkerSize',0.1);
h1=plot(glProfs.lon,glProfs.lat,'*m'); 
h2=plot(glProfs.argo.LONGITUDE,glProfs.argo.LATITUDE,'s','MarkerFaceColor','b','MarkerEdgeColor','k')
formatplot
legend([h1 h2],{'glider','ARGO'},'Location','SE')
title('(d)')

subtitle(t,'Sunfish and Argo Co-located Profiles (30 km, 5 days)')
t.TileSpacing='compact';
save_figure(gcf,['./plots/sunfish_argo_profiles'],[7.5 4.5],'.png','300')


%%
figure()
t = tiledlayout(2,2);


nexttile; hold on
for i = 1:length(glProfs.timeg)
    plot(glProfs.salinity(:,i)-glProfs.argo.PSAL_ADJUSTED(:,i),-glProfs.pressure(:,i),'Color',colMAP(i,:),'LineWidth',1.5);
end
ylim([-1000 0]); formatplot
xlabel('S^{glider} - S^{argo} / PSU')
title('(a)');

nexttile; hold on

for i = 1:length(glProfs.timeg)
    plot(glProfs.temperature(:,i)-glProfs.argo.TEMP_ADJUSTED(:,i),-glProfs.pressure(:,i),'Color',colMAP(i,:),'LineWidth',1.5);
end
ylim([-1000 0]); formatplot
xlabel('T^{glider} - T^{argo} / C')
formatplot
title('(b)')

nexttile; hold on

for i = 1:length(glProfs.timeg)
    plot(glProfs.oxygen_adjusted_umolkg(:,i)-glProfs.argo.DOX2_ADJUSTED(:,i),-glProfs.pressure(:,i),'Color',colMAP(i,:),'LineWidth',1.5);
end
ylim([-1000 0]); formatplot
xlabel('O^{glider}_2 - O^{argo}_2 / \mumol kg^{-1}')
formatplot
title('(c)')

nexttile; hold on
topoplot([-51.5 -49.5 50.7   54.9],[],[0:-1/2:-5],'-k')
plot(dat.gridded.lon,dat.gridded.lat,'.','MarkerSize',0.1);
h1=plot(glProfs.lon,glProfs.lat,'*m'); 
h2=plot(glProfs.argo.LONGITUDE,glProfs.argo.LATITUDE,'s','MarkerFaceColor','b','MarkerEdgeColor','k')
formatplot
legend([h1 h2],{'glider','ARGO'},'Location','SE')
title('(d)')

subtitle(t,'Difference Sunfish and Argo Co-located Profiles (30 km, 5 days)')
t.TileSpacing='compact';
save_figure(gcf,['./plots/sunfish_argo_profiles_diff'],[7.5 4.5],'.png','300')
