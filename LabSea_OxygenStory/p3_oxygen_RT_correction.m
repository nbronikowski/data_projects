clear; clc; close all;

%% WORK IN PROGRESS... Prob. do a monthly

toolbox_path = './processing_tools/';
addpath(genpath(toolbox_path), '-end');


path_name = './mat_files/'
var_name  = 'pearldiver_data'

load(fullfile(path_name,[var_name,'_oxy.mat']))

% annyoing to say pearldiver all the time
data = pearldiver;
clear pearldiver;

% select a few profiles
data.timeDateNum = datenum(data.time);
t1 = datenum(2020,02,10,0,0,0);
t2 = datenum(2020,02,20,0,0,0);
fn = fieldnames(data);
idx = find(data.timeDateNum>=t1 & data.timeDateNum<=t2);

for i = 1:length(fn)
    data.(fn{i}) = data.(fn{i})(idx); %#ok<*FNDSB> 
end

% Remove empty profiles
uidx = unique(data.profile_index);
uidx = uidx(mod(uidx,1)==0);
idx = (mod(data.profile_index,1)==0);
data.profile_index(~idx) =NaN;

for i = 1:length(uidx)
    idx = find(data.profile_index == uidx(i));
    els = data.profile_index(idx);
    if length(els(~isnan(els)))<10
        data.profile_index(idx)=NaN;
        data.profile_index(idx);
    end
end

data.oxygen_concentration_corr = NaN*data.raw_oxygen_concentration;

timeVec = data.timeDateNum;
timeVecSec = (timeVec-timeVec(1))*86400;

% Define your 1-second time grid
t_grid = (nanmin(timeVecSec):1:nanmax(timeVecSec))'; 

% Interpolate data onto the 1-second grid
idnan = ~isnan(data.oxygen_calphase_cleaned);
optode_phase_interp = interp1gap(timeVecSec(idnan), data.oxygen_calphase_cleaned(idnan), t_grid,40);
temperature_interp = interp1(timeVecSec(idnan), data.temperature(idnan), t_grid);
pressure_interp = interp1(timeVecSec(idnan), data.pressure(idnan), t_grid);
salinity_interp = interp1(timeVecSec(idnan), data.salinity(idnan), t_grid);
prof_idx_interp = interp1(timeVecSec(idnan), data.profile_index(idnan), t_grid,'nearest');
prof_dir_interp = interp1(timeVecSec(idnan), data.profile_direction(idnan), t_grid,'nearest');

% Define time scales
time_scales = 0:4:50;
tau_fac_fwd = 0;
tau_fac_bkwd = 0;

% Preallocate
virtual_foil_temps = zeros(length(temperature_interp), length(time_scales));

% Apply the exponential filter with different time constants
for i = 1:length(time_scales)
    virtual_foil_temps(:, i) = exponential_filter(temperature_interp,t_grid,time_scales(i),tau_fac_fwd);  % replace with your actual filter function
end

oxygen_concentrations = zeros(size(virtual_foil_temps));
foilcoef = [2.658710E-03,1.305085E-04,2.536294E-06,6.186505E-02,5.658259E-05,-1.514249E-02,1.192067E-03];
modeltype = 'uchidaAADI';
for i = 1:size(virtual_foil_temps, 2)
    oxygen_concentrations(:, i) = optcalcO2(virtual_foil_temps(:,i),optode_phase_interp, foilcoef, modeltype, salinity_interp, 1013.25, pressure_interp);
end

time_scales_reverse = 10:20:200;
oxygen_concentrations_filtered = zeros([size(oxygen_concentrations) length(time_scales_reverse)]);

for i = 1:length(time_scales_reverse)
    for j = 1:size(virtual_foil_temps, 2)
%         Apply your reverse exponential filter here
        oxygen_concentrations_filtered(:, j, i) = reverse_exponential_filter(oxygen_concentrations(:, j), t_grid, time_scales_reverse(i), tau_fac_bkwd);
    end
end

uprofs = unique(prof_idx_interp(~isnan(prof_idx_interp)));
RMSprof=NaN(size(virtual_foil_temps, 2),length(time_scales_reverse));

for i = 1:length(time_scales_reverse)
    for j = 1:size(virtual_foil_temps, 2)
        RMSprof(j,i) = nanmedian(calcUpDownRMS(oxygen_concentrations_filtered(:,j,i),t_grid,pressure_interp,prof_idx_interp,prof_dir_interp));
    end
end

time_scales_reverse_fine = time_scales_reverse(1):0.5:time_scales_reverse(end);
time_scales_fine = time_scales(1):0.5:time_scales(end);
[xq,yq]=meshgrid(time_scales_reverse_fine,time_scales_fine);
RMSprof_fine = interp2(time_scales_reverse,time_scales,RMSprof,xq,yq);
[~, minIndex] = nanmin(RMSprof_fine(:));
[idx_rev_time, idx_time_scale] = ind2sub(size(RMSprof_fine), minIndex);
rev_time =  time_scales_reverse_fine(idx_time_scale)
adv_time = time_scales_fine(idx_rev_time)


figure(); contourf(xq,yq,RMSprof_fine); cb=colorbar; hold
xlabel('\tau_{O_2} (s)')
ylabel('\tau_T (s)')
plot(rev_time,adv_time,'+r','MarkerSize',20,'LineWidth',5)
ylabel(cb,'RMSE')
colormap(gca,cmocean('-solar',15));
% caxis([0 max(RMSprof(:))])
title('ADI 4831 SN 124 Slow Foil Response Time Constants vs Profile RMS')
save_figure(gcf,'pearldiver_oxy_response',[6 5],'.png','300')


% Apply to oxygen data
foil_temps = exponential_filter(temperature_interp,t_grid,adv_time,tau_fac_fwd);  % replace with your actual filter function

foilcoef = [2.658710E-03,1.305085E-04,2.536294E-06,6.186505E-02,5.658259E-05,-1.514249E-02,1.192067E-03];
modeltype = 'uchidaAADI';
oxygen_concentrations = optcalcO2(foil_temps,optode_phase_interp, foilcoef, modeltype, salinity_interp, 1013.25, pressure_interp);
oxygen_concentrations_filtered = reverse_exponential_filter(oxygen_concentrations,t_grid,rev_time,tau_fac_bkwd);  % replace with your actual reverse filter function

[b, a] = butter(3, 0.05); % 3rd order butterworth filter like in GEOMAR toolbox
oxygen_concentrations_raw   = optcalcO2(temperature_interp,optode_phase_interp, foilcoef, modeltype, salinity_interp, 1013.25, pressure_interp);
oxygen_concentrations_filt  = safe_filtfilt(b,a,oxygen_concentrations_filtered);

figure(); hold on
plot(t_grid,oxygen_concentrations_raw,'b')
plot(t_grid,oxygen_concentrations_filt,'r')

%     data.oxygen_concentration_corr(seg_idx) = interp1gap(t_grid,oxygen_concentrations_filt,timeVecSec,30);
%     
%     %% Gridding and interpolation
%     pg = 0:1:ceil(nanmax(pressure_interp)); 
%     [~,~,pg_time] = pgrid_columns(prof_idx_interp,pressure_interp,t_grid,pg);
%     [~,~,pg_o2conc] = pgrid_columns(prof_idx_interp,pressure_interp,oxygen_concentrations_filt,pg);
%     
%     pg_o2conc_filled = inpaint_nans(pg_o2conc,3);  % Interpolate using piecewise cubic Hermite interpolating polynomial
%     pg_o2conc_filt = medfilt2(pg_o2conc_filled,[1 2],'symmetric');
%     pg_o2conc_filt(isnan(pg_o2conc)) = NaN;
% 
% %     figure(); hold on
% %     h1=plot(optode_conc_interp,-pressure_interp,'.','Color',[.6 .6 .6]);
% %     h2=plot(pg_o2conc_filt,-pg,'-k');
% %     legend({'raw sci\_ox4 (salinity compensated)','response-time corrected'},...
% %         'location','SW');
% %     formatplot
% %     ylim([-max(pg) 0])
% 
% %     ylabel('Pressure (dbar)')
% %     xlabel('Oxygen Concentration (\mumol L^{-1})')
% %     title('Aanderaa Optode 4831 w. Fast Foil Correction')
% %     save_figure(gcf,'oxy_correction',[6 5],'.png','300')
% % 
% %     pg_o2conc_filt_ds=spike_test(pg_o2conc_filt,100,3,0);
% %     figure()
% %     subplot(121)
% %     imagescn(pg_o2conc)
% %     
% %     subplot(122)
% %     imagescn(pg_o2conc_filt_ds)
% 
%     
%     gridded_oxy  = pg_o2conc_filt(:);
%     time_oxy = pg_time(:);
%     [time_oxy,idor] = sort(time_oxy);
%     gridded_oxy = gridded_oxy(idor);
%     [time_oxy,idu]=unique(time_oxy);
%     gridded_oxy = gridded_oxy(idu);
%     data.oxygen_concentration_gridded_corr(seg_idx) =interp1gap(time_oxy,gridded_oxy,timeVecSec,30);
% %     
% %     


% data.oxygen_solubility = O2solubility(data.temperature,data.salinity_corrected);
% data.oxygen_saturation = data.oxygen_concentration_gridded_corr./data.oxygen_solubility*100;
% save('glider_data_oxy_processed.mat','data');







