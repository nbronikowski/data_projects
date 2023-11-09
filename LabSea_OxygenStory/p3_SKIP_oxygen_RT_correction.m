clear; clc; close all;

% Update, Nov 9:  
% -------------------------------------------------------------------------
% Suggesting we skip this method at least for Pearldiver
% Reason is that gradients in temperature are low and we are imparting 
% significant noise to the glider oxygen data without necessarily gaining 
% confidence in data.
%
% Conclusion:
% Suggesting we just compare with mooring data in density and do an in-situ 
% calibration in density.


toolbox_path = './processing_tools/';
addpath(genpath(toolbox_path), '-end');


path_name = './mat_files/'
var_name  = 'pearldiver_data'

load(fullfile(path_name,[var_name,'_oxy.mat']))

% annyoing to say pearldiver all the time
data = pearldiver;
clear pearldiver;

%% Remove empty profiles caused by init of house elf at surface
uidx = unique(data.profile_index);
uidx = uidx(mod(uidx,1)==0);
idx = (mod(data.profile_index,1)==0);
data.profile_index(~idx) =NaN;
for i = 1:length(uidx)
    idx = find(data.profile_index == uidx(i));
    els = data.profile_index(idx);
    if length(els(~isnan(els)))<10
        data.profile_index(idx)=NaN;
    end
end

%% Select a Time Subset
data.timeDateNum = datenum(data.time);
% t1 = datenum(2020,05,10,0,0,0);
% t2 = datenum(2020,05,15,0,0,0);
% fn = fieldnames(data);
% idx = find(data.timeDateNum>=t1 & data.timeDateNum<=t2);
% for i = 1:length(fn)
%     data.(fn{i}) = data.(fn{i})(idx); %#ok<*FNDSB> 
% end

%% Loop every day of year
doy=doy(data.timeDateNum);
doy_sf = scale_var(doy,1);
udd = unique(doy_sf);

data.oxygen_concentration_corr = NaN*data.raw_oxygen_concentration;

for ii = 1:length(udd)

    % get doy subset (data for glider that day)
    idx_time_seg = doy_sf == udd(ii);

    idnan = find(~isnan(data.oxygen_calphase_cleaned(idx_time_seg)));
    if length(idnan)>1000 % should be real data in oxygen to do this correction
        timeVec = data.timeDateNum(idx_time_seg);
        optode_phase = data.oxygen_calphase_cleaned(idx_time_seg);
        temperature  = data.temperature(idx_time_seg);
        pressure     = data.pressure(idx_time_seg);
        salinity     = data.salinity(idx_time_seg);
        profile_index= data.profile_index(idx_time_seg);
        profile_direction= data.profile_direction(idx_time_seg);
    
        % Interpolate data onto the 1-second grid but only use data w. no nan's
        % in oxygen phase
        
        % Define your 1-second time grid
        timeVecSec = (timeVec-timeVec(1))*86400;
        t_grid = (nanmin(timeVecSec):1:nanmax(timeVecSec))'; 
    
        optode_phase_interp = interp1gap(timeVecSec(idnan), optode_phase(idnan), t_grid,40);
        temperature_interp = interp1(timeVecSec(idnan), temperature(idnan), t_grid);
        pressure_interp = interp1(timeVecSec(idnan), pressure(idnan), t_grid);
        salinity_interp = interp1(timeVecSec(idnan), salinity(idnan), t_grid);
        prof_idx_interp = interp1(timeVecSec(idnan), profile_index(idnan), t_grid,'nearest');
        prof_dir_interp = interp1(timeVecSec(idnan), profile_direction(idnan), t_grid,'nearest');
        
        % Define time scales
        time_scales = 0:4:50;
        time_scales_reverse = 0:10:200;
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
        rev_time(ii) =  time_scales_reverse_fine(idx_time_scale);
        adv_time(ii) = time_scales_fine(idx_rev_time);
        
        %% Plot RT Fitting Result
        % figure(); contourf(xq,yq,RMSprof_fine); cb=colorbar; hold
        % xlabel('\tau_{O_2} (s)')
        % ylabel('\tau_T (s)')
        % plot(rev_time,adv_time,'+r','MarkerSize',20,'LineWidth',5)
        % ylabel(cb,'RMSE')
        % colormap(gca,cmocean('-solar',15));
        % % caxis([0 max(RMSprof(:))])
        % title('ADI 4831 SN 124 Slow Foil Response Time Constants vs Profile RMS')
        % save_figure(gcf,'./plots/pearldiver_oxy_response',[6 5],'.png','300')
        
        % Apply to oxygen data
        foil_temps = exponential_filter(temperature_interp,t_grid,adv_time(ii),tau_fac_fwd);  % replace with your actual filter function
        oxygen_concentrations = optcalcO2(foil_temps,optode_phase_interp, foilcoef, modeltype, salinity_interp, 1013.25, pressure_interp);
        oxygen_concentrations_filtered = reverse_exponential_filter(oxygen_concentrations,t_grid,rev_time(ii),tau_fac_bkwd);  % replace with your actual reverse filter function
        
        [b, a] = butter(3, 0.05); % 3rd order butterworth filter like in GEOMAR toolbox
        oxygen_concentrations_raw   = optcalcO2(temperature_interp,optode_phase_interp, foilcoef, modeltype, salinity_interp, 1013.25, pressure_interp);
        oxygen_concentrations_filt  = safe_filtfilt(b,a,oxygen_concentrations_filtered);
         
        %% Plot RT Corrected Data
        %figure(); hold on
        %plot(t_grid,oxygen_concentrations_raw,'b')
        %plot(t_grid,oxygen_concentrations_filt,'r')
        disp(num2str(ii/length(udd)*100))

        data.oxygen_concentration_corr(idx_time_seg)= interp1gap(t_grid,oxygen_concentrations_filt,timeVecSec,40);
    end
end
data.oxygen_concentration = optcalcO2(data.temperature,...
    data.oxygen_calphase_cleaned,foilcoef,modeltype,data.salinity,1013,data.pressure);
save('oxy_tau_vals_pearldiver.mat','rev_time','adv_time')


