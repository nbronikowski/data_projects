clear; clc; close all;

path_name = './mat_files/'
var_name  = 'pearldiver_data'
load(fullfile(path_name,[var_name,'_oxy_qc.mat']))
dat = pearldiver; 
dat.gridded.pressure_grid=ones(size(dat.gridded.pressure)).*dat.gridded.pressure_grid';

dat.gridded.rho = sw_dens(dat.gridded.salinity,dat.gridded.temperature,dat.gridded.pressure);
dat.gridded.oxygen_adjusted_umolkg = dat.gridded.oxygen_adjusted./(dat.gridded.rho/1000);
dat.gridded.oxygen_raw_umolkg = dat.gridded.oxygen_raw./(dat.gridded.rho/1000);

%% Load Argo Data
load('./ARGO/argo_oxygen.mat')

% Co-locate glider and Argo measurements and store result
idx = find(argo.TIME>=dat.gridded.timeg(1) & argo.TIME<=dat.gridded.timeg(end));

GL_VAR_LIST = fieldnames(dat.gridded);
AR_VAR_LIST = fieldnames(argo);
argo.pgrid = ones(size(argo.PRES)).*argo.pgrid;
m= 1;
for i = 1:length(argo.TIME(idx))
    dkm=lldistkm([dat.gridded.lat(:),dat.gridded.lon(:)],...
        [argo.LATITUDE(idx(i)),argo.LONGITUDE(idx(i))]);
    dtime = abs(dat.gridded.timeg - argo.TIME(idx(i)));
    dkm=dkm(:); dtime = dtime(:);
    dkm(dkm>30)=NaN; 
    dtime(dtime>5)=NaN;
    [val,glIDX] = min(dkm(:)+dtime(:),[],'omitnan');
    if ~isempty(glIDX) && val<30

        for j = 1:length(GL_VAR_LIST)
            glProfs.(GL_VAR_LIST{j})(:,m)=dat.gridded.(GL_VAR_LIST{j})(:,glIDX);
        end
        for k = 1:length(AR_VAR_LIST)
            glProfs.argo.(AR_VAR_LIST{k})(:,m) = argo.(AR_VAR_LIST{k})(:,idx(i));
        end
        m = m+1;
    end
end
tARGO = ones(size(glProfs.argo.TEMP)).*glProfs.argo.TIME;



%% PLOT RESULTING STRUCTURE

figure()
t = tiledlayout(2,2);


nexttile; hold on
colMAP = seminfhaxby(length(glProfs.timeg)+2);
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
topoplot([-52.8  -50.8 56.   57.4],[],[0:-1/2:-5],'-k')
h1=plot(glProfs.lon,glProfs.lat,'*m'); 
h2=plot(glProfs.argo.LONGITUDE,glProfs.argo.LATITUDE,'s','MarkerFaceColor','b','MarkerEdgeColor','k')
formatplot
legend([h1 h2],{'glider','ARGO'},'Location','SE')
title('(d)')

subtitle(t,'Difference Pearldiver and Argo Co-located Profiles (30 km, 5 days)')
t.TileSpacing='compact';
save_figure(gcf,['./plots/pearldiver_argo_profiles'],[7.5 4.5],'.png','300')


%%
figure()
t = tiledlayout(2,2);


nexttile; hold on
colMAP = seminfhaxby(length(glProfs.timeg)+2);
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
topoplot([-52.8  -50.8 56.   57.4],[],[0:-1/2:-5],'-k')
h1=plot(glProfs.lon,glProfs.lat,'*m'); 
h2=plot(glProfs.argo.LONGITUDE,glProfs.argo.LATITUDE,'s','MarkerFaceColor','b','MarkerEdgeColor','k')
formatplot
legend([h1 h2],{'glider','ARGO'},'Location','SE')
title('(d)')

subtitle(t,'Difference Pearldiver and Argo Co-located Profiles (30 km, 5 days)')
t.TileSpacing='compact';
save_figure(gcf,['./plots/pearldiver_argo_profiles_diff'],[7.5 4.5],'.png','300')

