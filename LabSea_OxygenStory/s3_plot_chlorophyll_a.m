clear; clc; close all;

path_name = './mat_files/'
var_name  = 'sunfish_data'
load(fullfile(path_name,[var_name,'_clean.mat']))
dat = sunfish;  clear sunfish

addpath('./processing_tools/');

% % Dark Counts Calc. in-situ
% data_range = dat.pressure>350 & dat.pressure<900 & ...
%     dat.time>datetime(2022,01,20) & dat.time<datetime(2022,04,01);
% 
% %% Backscatter 700nm 
% 
% eco_bb_dark_counts = 47; % factory value
% eco_bb_scale_factor = 1.913E-06; % m-1 sr-1 / counts
% % eco_bb_dark_counts = 110; %min_75th_percentile(dat.bb700_sig(data_range))
% 
% lambda = 700; % 700 nm backscatter 
% theta = 117; % angle of instrument
% VBSC  = eco_bb_scale_factor*(dat.bb700_sig-eco_bb_dark_counts);  
% 
% % 1. Method Zhang et al., 2009 Volume Scattering of Seawater
% beta_sw=betasw(lambda, dat.temperature, theta, dat.salinity_cor);
% 
% % Particle backscatter coefficient, BBP 
% Xp = 1.1; %recommended by WETLabs 2013
% dat.bbp700 = 2*pi*(VBSC - beta_sw) * Xp;
% 
% %% Chlorophyll
% eco_chl_dark_counts = 45; % factory value
% eco_chl_scale_factor = 0.0073; % micro-g L-1
% 
% % Investigate chl_dark_counts from values when readings where dark
% eco_chl_dark_counts = min_75th_percentile(dat.chlor_sig(data_range))
% 
% dat.chlor  = eco_chl_scale_factor*(dat.chlor_sig-eco_chl_dark_counts);  
% 
% %% Cdom
% eco_cdom_dark_counts = 50;
% eco_cdom_scale_factor = 0.0904;
% eco_cdom_dark_counts = min_75th_percentile(dat.cdom_sig(data_range))
% 
% dat.cdom  = eco_cdom_scale_factor*(dat.cdom_sig-eco_cdom_dark_counts);  
% 
% id = abs(dat.roll)>5 | dat.chlor<0 | dat.bbp700<0 | dat.cdom<0;
% 
% dat.chlor(id)=NaN;
% dat.bbp700(id)=NaN;
% dat.cdom(id)=NaN;
% 
% % next grid chlorophyll, cdom and bb700 
% pg = 0:1:ceil(max(dat.pressure,[],'omitnan')); 
% 
% [~,~,conduc] = pgrid_columns(dat.profile_index,dat.pressure,dat.conductivity,pg);
% [~,~,chlor] = pgrid_columns(dat.profile_index,dat.pressure,dat.chlor,pg);
% [~,~,cdom] = pgrid_columns(dat.profile_index,dat.pressure,dat.cdom,pg);
% [~,~,bbp700] = pgrid_columns(dat.profile_index,dat.pressure,dat.bbp700,pg);
% 
% [~,column_idx]=deleteAlmostEmptyColumns(conduc,pg);
% dat.gridded.chlor = chlor(:,column_idx);
% dat.gridded.cdom  = cdom(:,column_idx);
% dat.gridded.bbp700 = bbp700(:,column_idx);
% 
% 
% dat.gridded.timeg = mean(dat.gridded.time,1,'omitnan');
% dat.gridded.lon = interp1(datenum(dat.time),dat.longitude,dat.gridded.timeg);
% dat.gridded.lat = interp1(datenum(dat.time),dat.latitude,dat.gridded.timeg);


%% Remove Empty Columns
[pg_chlor,column_idxCHLOR]=deleteAlmostEmptyColumns(dat.gridded.chlor,...
    dat.gridded.pressure_grid,20);
pg_bbp700 = dat.gridded.bbp700(:,column_idxCHLOR);
pg_time = dat.gridded.timeg(1,column_idxCHLOR);
pg_lat  = dat.gridded.lat(1,column_idxCHLOR);

pg_bbp700_filt=medfilt2(pg_bbp700,[4 4],'symmetric');
dx = (pg_bbp700-pg_bbp700_filt);
pg_chlor(abs(dx)>0.002 )=NaN;
pg_bbp700(abs(dx)>0.002 )=NaN;

[a_zenith,a_zenith_rad] = solar_zenith(pg_time,pg_lat);
a_zenith = a_zenith(:);
pg_time = pg_time(:);
pg = dat.gridded.pressure_grid;

pg = pg(:);

%% Quenching Correction
thomalla_correction = quenching_correction(a_zenith,pg_time,pg_chlor,pg_bbp700);

chla_corr = [];
chla_night = [];

for i = 1:length(thomalla_correction)
    temp_corrected =  thomalla_correction(i).thomallacorrect;
    temp_night     = thomalla_correction(i).chl_night;
    if isempty(temp_corrected)
        temp_corrected = thomalla_correction(i).chl_day;
    end
    
    temp = [temp_corrected,temp_night];
    chla_corr = [chla_corr,temp];
    temp = [temp_corrected*NaN,thomalla_correction(i).chl_night];
    chla_night = [chla_night,temp];
end


figure()
subplot(311)
imagescn(chla_night);
set(gca,'ydir','reverse'); caxis([0 2]); colorbar
ylim([0 400])

subplot(312)
imagescn(chla_corr);
set(gca,'ydir','reverse'); caxis([0 2]); colorbar
ylim([0 400])

subplot(313)
imagescn(pg_bbp700);
set(gca,'ydir','reverse'); colorbar
caxis([0 0.01])
ylim([0 400])

%% Integrate Chlorophyll Vertically
intg_chla = nan(length(pg_time), 1);  % Preallocate the result vector
intg_corr_chla = nan(length(pg_time), 1);  % Preallocate the result vector
intg_night_chla = nan(length(pg_time),1);
z_thrs = 30;
for i = 1:length(pg_time)
    var_nan_idx = ~isnan(pg_chlor(:,i));
    if abs(max(pg(var_nan_idx))-min(pg(var_nan_idx)))>z_thrs

        idx_start = find(pg(var_nan_idx)>0,1,'first'); % from 0
        idx_end = find(pg(var_nan_idx)<100,1,'last');
        
        temp_z = pg(idx_start:idx_end);
        temp_y = pg_chlor(idx_start:idx_end,i);
        temp_y2 = chla_corr(idx_start:idx_end,i);
        temp_y3 = chla_night(idx_start:idx_end,i);
        
        % interp 
        idnan = find(~isnan(temp_y));
        idnan2 = find(~isnan(temp_y2));
        temp_y = interp1(temp_z(idnan),temp_y(idnan),temp_z,'linear','extrap');

        if isempty(idnan2)
            temp_y2 = temp_z*NaN;
        else
            temp_y2 = interp1(temp_z(idnan2),temp_y2(idnan2),temp_z,'linear','extrap');
        end


        % Flip vectors if depth is in descending order
        if temp_z(1) > temp_z(end)
            temp_z = flip(temp_z);
            temp_y = flip(temp_y);
            temp_y2 = flip(temp_y2);
            temp_y3 = flip(temp_y3);
        end
        
        % Integrate
        intg_chla(i) = trapz(temp_z, temp_y);
        intg_corr_chla(i) = trapz(temp_z, temp_y2);

        if length(temp_y3(~isnan(temp_y3)))>4
            idnan3 = find(~isnan(temp_y3));
            temp_y3 = interp1(temp_z(idnan3),temp_y3(idnan3),temp_z,'linear','extrap');
            intg_night_chla(i) = trapz(temp_z,temp_y3);
        else
            intg_night_chla(i) = NaN;
        end
        

    end

end

%% Plot Result

figure()
t=tiledlayout(3,1);
t1 = datenum(2022,1,15);
t2 = datenum(2022,5,20);

nexttile
imagescn(pg_time,-pg,pg_chlor)
xlim([t1 t2])
ylim([-400 0])
ylabel('depth / m')
cb=colorbar; colormap(gca,cmocean('algae'))
caxis([0 2])
title('a) Chlorophyll-a Concentrations')
formatplot
ylabel(cb,'Chl_a / mg m^{-3}')
datetick('x','dd-mmm-yy','keeplimits')

nexttile
imagescn(pg_time,-pg,chla_corr)
xlim([t1 t2])
ylim([-400 0])
ylabel('depth / m')
cb2=colorbar; colormap(gca,cmocean('algae'))
caxis([0 2])
title('b) Quenching Corrected Chlorophyll-a Concentrations')
formatplot
ylabel(cb2,'Chl_a / mg m^{-3}')
datetick('x','dd-mmm-yy','keeplimits')

nexttile; hold on
plot(pg_time,intg_corr_chla,'-b')
plot(pg_time,intg_chla,'-r')
legend('Chl_a Thomalla (2017)','Chl_a (raw)','Location','best')
xlim([t1 t2])
ylabel('Chl^{[0<z<200m]}_a / mgm^{-3}')
title('c) Vertically Integrated Chloropyhll (0-200 m)')
datetick('x','dd-mmm-yy','keeplimits')
formatplot


t.TileSpacing='compact';
t.Padding='compact';
% subtitle(t,'Sunfish Glider Chlorophyll-a Concentrations')

save_figure(gcf,['./plots/sunfish_chlorophyll_a'],[7.5 6],'.png','300')
