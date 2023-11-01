clear; close all; clc
load glider_data_CTD_corrected.mat

%% Backscatter Correction
% Compute Volume Backscatter
eco_bb_dark_counts = 49;
eco_bb_scale_factor = 1.693E-06; % m-1 sr-1 / counts
lambda = 700; % 700 nm backscatter 
theta = 117; % angle of instrument
VBSC  = eco_bb_scale_factor*(gdat.bb700_sig-eco_bb_dark_counts);  

% 1. Method Zhang et al., 2009
% Volume Scattering of Seawater
beta_sw=betasw(lambda, gdat.temp, theta, gdat.salt_corr);

% Particle backscatter coefficient, BBP 
Xp = 1.1; %recommended by WETLabs 2013
gdat.bbp700 = 2*pi*(VBSC - beta_sw) * Xp;

% 2. Method from Wetlabs
% Volume Scattering of particles adjusted to seawater
% delta = 0.09; % from Wetlabs
% beta_w = 1.38*(lambda/500)^-4.32*(1+0.3.*gdat.salt_corr/ 37)*...
%     1e-4*(1+cosd(theta)^2*(1-delta)/(1+delta));

% % Scattering data correction - Ignoring this but found in the Manual!
% a = 1;
% VBSC_corr = VBSC*exp(0.0391*a);
% 
% % compute volume scattering corrected for seawater particles
% beta_p = VBSC - beta_w;

% % Backscattering coefficient from single measurement w/out effect for
% seawater
% bbp700 = 2*pi*Xp*beta_p;

%% Chlorophyll
eco_chl_dark_counts = 49; % factory value

% Investigate chl_dark_counts
gdat.water_depth(gdat.water_depth<10)=NaN;
idnan = ~isnan(gdat.water_depth);
gdat.water_depth = interp1gap(gdat.timeDateNum(idnan),gdat.water_depth(idnan),gdat.timeDateNum,2000/86400);
eco_chl_dark_counts = nanmin(gdat.chlor_sig(gdat.pres<170 & gdat.pres>150 & abs(gdat.depth-gdat.water_depth)>30));

eco_chl_scale_factor = 0.0072; % micro-g L-1


gdat.chla  = eco_chl_scale_factor*(gdat.chlor_sig-eco_chl_dark_counts);  


gdato = gdat;

% select a few profiles
t1 = datenum(2023,08,20,0,0,0);
t2 = datenum(2023,08,30,0,0,0);
fn = fieldnames(gdat);
idx = find(gdat.timeDateNum>=t1 & gdat.timeDateNum<=t2);
for i = 1:length(fn)
    gdat.(fn{i}) = gdat.(fn{i})(idx); %#ok<*FNDSB> 
end

% gridding routine based on u_prof_idx
pg = 0:1:200; 

% calc sigma0
gdat.sigma_t0=sw_dens0(gdat.salt_corr,gdat.temp)-1000; % sigma-t wrp to 0 dbar

% pgrid in coloums
[~,~,pg_chla] = pgrid_columns(gdat.prof_idx,gdat.pres,gdat.chla,pg);
[~,~,pg_salt] = pgrid_columns(gdat.prof_idx,gdat.pres,gdat.salt_corr,pg);
[~,~,pg_sigma_t0] = pgrid_columns(gdat.prof_idx,gdat.pres,gdat.sigma_t0,pg);
[~,~,pg_lat] = pgrid_columns(gdat.prof_idx,gdat.pres,gdat.lat,pg);
[~,~,pg_topo_z] = pgrid_columns(gdat.prof_idx,gdat.pres,gdat.water_depth,pg);
[~,~,pg_lon] = pgrid_columns(gdat.prof_idx,gdat.pres,gdat.lon,pg);
[~,~,pg_temp] = pgrid_columns(gdat.prof_idx,gdat.pres,gdat.temp,pg);
[~,~,pg_bbp700] = pgrid_columns(gdat.prof_idx,gdat.pres,gdat.bbp700,pg);
[~,~,pg_time] = pgrid_columns(gdat.prof_idx,gdat.pres,gdat.timeDateNum,pg);

[pg_chla,idx]=deleteAlmostEmptyColumns(pg_chla,pg);
pg_temp = pg_temp(:,idx);
pg_salt = pg_salt(:,idx);
pg_bbp700 = pg_bbp700(:,idx);
pg_sigma_t0 = pg_sigma_t0(:,idx);
pg_time = pg_time(:,idx);
pg_lat   = pg_lat(:,idx);
pg_lon   = pg_lon(:,idx);
pg_topo_z   = pg_topo_z(:,idx);


timeg  = nanmean(pg_time,1);
long  = nanmean(pg_lon,1);
latg  = nanmean(pg_lat,1);
z_wat = nanmean(pg_topo_z,1);
D1 = 1:length(timeg); idn = ~isnan(timeg);
timeg  = interp1(D1(idn),timeg(idn),D1);
[x,y] = meshgrid(timeg,pg);


%% Check quenching correction
[a_zenith,a_zenith_rad] = solar_zenith(timeg,latg);

thomalla_correction = quenching_correction(a_zenith,timeg,pg_chla,pg_bbp700);

chla_corr = [];
chla_night = [];
for i = 1:length(thomalla_correction)
    temp = [thomalla_correction(i).thomallacorrect,thomalla_correction(i).chl_night];
    chla_corr = [chla_corr,temp];
    temp = [thomalla_correction(i).thomallacorrect*NaN,thomalla_correction(i).chl_night];
    chla_night = [chla_night,temp];
end

%% Compute mixed layer abs() density diff approach 0.05 kg m3
pmld_dens = mixed_layer(pg_salt,pg_temp,pg,pg_sigma_t0,0.05);
% %% Holte & Talley method
% % for i = 1:length(timeg)
% %     pmld_opt_ht(:,i) = mld(pg(:),pg_temp(:,i),pg_salt(:,i));
% % end

plot_name = 'thomalla_correction';
% p2_analysis_plots
%% It looks like quenching correction is not really helping due to cloud cover...
% so we skip it in further analysis


%% Vertically integrated chlorophyll 
intg_chla = nan(length(timeg), 1);  % Preallocate the result vector
intg_corr_chla = nan(length(timeg), 1);  % Preallocate the result vector
intg_night_chla = nan(length(timeg),1);
z_thrs = 50;
for i = 1:length(timeg)
    var_nan_idx = ~isnan(pg_chla(:,i));
    if abs(max(pg(var_nan_idx))-min(pg(var_nan_idx)))>z_thrs

        idx_start = find(pg(var_nan_idx)>pmld_dens(i), 1 ); % from below mld to z_thrs m + mld
        idx_end = find(pg(var_nan_idx)<z_thrs+pmld_dens(i), 1, 'last' ); % z_thrs m stop from mld
        
        temp_z = pg(idx_start:idx_end);
        temp_y = pg_chla(idx_start:idx_end,i);
        temp_y2 = chla_corr(idx_start:idx_end,i);
        temp_y3 = chla_night(idx_start:idx_end,i);
        

        % interp 
        idnan = find(~isnan(temp_y));
        idnan2 = find(~isnan(temp_y2));
        temp_y = interp1(temp_z(idnan),temp_y(idnan),temp_z,'linear','extrap');
        temp_y2 = interp1(temp_z(idnan2),temp_y2(idnan2),temp_z,'linear','extrap');
        
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

%% Plot Integrated Chlorophyll
plot_name = 'plot_integrated_chla';
% p2_analysis_plots

%% Plot CMEMS Multi Product 
plot_name = 'plot_CMEMS_chla';
% p2_analysis_plots

%% Plot Glider Section
plot_name = 'event_section';







