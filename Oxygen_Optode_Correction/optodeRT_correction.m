clear; clc; close all; clear all
load glider_data_processed.mat

data.oxygen_concentration_corr = NaN*data.oxygen_concentration;
data.oxygen_concentration_gridded_corr = NaN*data.oxygen_concentration;

idSeg = [1;find([0;diff(data.time)]>0.5);length(data.time)];
for kk = 1:length(idSeg)-1

    seg_idx = (idSeg(kk):idSeg(kk+1)-1)'; % gives the section index
    
    timeVec = data.time(seg_idx);
    timeVecSec = (timeVec-timeVec(1))*86400;

    % Define your 1-second time grid
    t_grid = (nanmin(timeVecSec):1:nanmax(timeVecSec))'; 
    
    % Interpolate data onto the 1-second grid
    optode_phase_interp = interp1gap(timeVecSec, data.oxygen_calphase(seg_idx), t_grid,10);
    optode_conc_interp = interp1gap(timeVecSec, data.oxygen_concentration(seg_idx), t_grid,10);
    optode_temp_interp = interp1gap(timeVecSec, data.oxygen_sensor_temperature(seg_idx), t_grid,10);
    temperature_interp= interp1gap(timeVecSec, data.temperature(seg_idx), t_grid,10);
    pressure_interp = interp1gap(timeVecSec, data.pressure(seg_idx), t_grid,10);
    salinity_interp = interp1gap(timeVecSec, data.salinity_corrected(seg_idx), t_grid,10);
    

    prof_idx_interp = interp1(timeVecSec, data.profile_index(seg_idx), t_grid,'nearest');
    prof_dir_interp = interp1(timeVecSec, data.profile_direction(seg_idx), t_grid,'nearest');


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
    foilcoef = [3.540882E-03	1.763056E-04	4.367122E-06	2.489593E+02	8.761237E-01	-5.441802E+01	4.447275E+00];
    modeltype = 'uchidaAADI';
    for i = 1:size(virtual_foil_temps, 2)
        oxygen_concentrations(:, i) = optcalcO2(virtual_foil_temps(:,i),optode_phase_interp, foilcoef, modeltype, salinity_interp, 1013.25, pressure_interp);
    end


    time_scales_reverse = 10:20:200;
    oxygen_concentrations_filtered = zeros([size(oxygen_concentrations) length(time_scales_reverse)]);
    
    for i = 1:length(time_scales_reverse)
        for j = 1:size(virtual_foil_temps, 2)
            % Apply your reverse exponential filter here
            oxygen_concentrations_filtered(:, j, i) = reverse_exponential_filter(oxygen_concentrations(:, j), t_grid, time_scales_reverse(i), tau_fac_bkwd);
        end
    end

    uprofs = unique(prof_idx_interp(~isnan(prof_idx_interp)));
    RMSprof=NaN(size(virtual_foil_temps, 2),length(time_scales_reverse));

    for i = 1:length(time_scales_reverse)
        for j = 1:size(virtual_foil_temps, 2)
            RMSprof(j,i) = nanmedian(calcRMS(oxygen_concentrations_filtered(:,j,i),t_grid,pressure_interp,prof_idx_interp,prof_dir_interp));
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


%     figure(); contourf(xq,yq,RMSprof_fine,'LevelStep',4); cb=colorbar; hold
%     xlabel('\tau_{O_2} (s)')
%     ylabel('\tau_T (s)')
%     plot(rev_time,adv_time,'+r','MarkerSize',20,'LineWidth',5)
%     ylabel(cb,'RMSE')
%     colormap(gca,cmocean('-solar',15));
%     caxis([0 max(RMSprof(:))])
%     title('ADI 4831 Fast Foil Response Time Constants vs Profile RMS')
%     save_figure(gcf,'oxy_response',[6 5],'.png','300')

        



    %% Apply to oxygen data
    foil_temps = exponential_filter(temperature_interp,t_grid,adv_time,tau_fac_fwd);  % replace with your actual filter function

    foilcoef = [3.540882E-03	1.763056E-04	4.367122E-06	2.489593E+02	8.761237E-01	-5.441802E+01	4.447275E+00];
    modeltype = 'uchidaAADI';
    oxygen_concentrations = optcalcO2(foil_temps,optode_phase_interp, foilcoef, modeltype, salinity_interp, 1013.25, pressure_interp);
    oxygen_concentrations_filtered = reverse_exponential_filter(oxygen_concentrations,t_grid,rev_time,tau_fac_bkwd);  % replace with your actual reverse filter function
   
    [b, a] = butter(3, 0.05); % 3rd order butterworth filter like in GEOMAR toolbox

    % plot and compare different models
    %     figure(); hold on
    %     plot(t_grid,oxygen_concentrations_filtered,'b')
    %     plot(t_grid,movmedian(oxygen_concentrations_filtered,20,1,'omitnan'),'LineWidth',2,'Color','k')
    %     plot(t_grid,safe_filtfilt(b,a,oxygen_concentrations_filtered),'LineWidth',2,'Color','r')
    %     plot(t_grid,sgolayfilt(oxygen_concentrations_filtered,3,41),'LineWidth',2,'Color','m')
    %     plot(t_grid,sgolayfilt(oxygen_concentrations_filtered,3,41),'LineWidth',2,'Color','m')
    
    oxygen_concentrations_filt  = safe_filtfilt(b,a,oxygen_concentrations_filtered);
    data.oxygen_concentration_corr(seg_idx) = interp1gap(t_grid,oxygen_concentrations_filt,timeVecSec,30);
    
    %% Gridding and interpolation
    pg = 0:1:ceil(nanmax(pressure_interp)); 
    [~,~,pg_time] = pgrid_columns(prof_idx_interp,pressure_interp,t_grid,pg);
    [~,~,pg_o2conc] = pgrid_columns(prof_idx_interp,pressure_interp,oxygen_concentrations_filt,pg);
    
    pg_o2conc_filled = inpaint_nans(pg_o2conc,3);  % Interpolate using piecewise cubic Hermite interpolating polynomial
    pg_o2conc_filt = medfilt2(pg_o2conc_filled,[1 2],'symmetric');
    pg_o2conc_filt(isnan(pg_o2conc)) = NaN;

%     figure(); hold on
%     h1=plot(optode_conc_interp,-pressure_interp,'.','Color',[.6 .6 .6]);
%     h2=plot(pg_o2conc_filt,-pg,'-k');
%     legend({'raw sci\_ox4 (salinity compensated)','response-time corrected'},...
%         'location','SW');
%     formatplot
%     ylim([-max(pg) 0])

%     ylabel('Pressure (dbar)')
%     xlabel('Oxygen Concentration (\mumol L^{-1})')
%     title('Aanderaa Optode 4831 w. Fast Foil Correction')
%     save_figure(gcf,'oxy_correction',[6 5],'.png','300')
% 
%     pg_o2conc_filt_ds=spike_test(pg_o2conc_filt,100,3,0);
%     figure()
%     subplot(121)
%     imagescn(pg_o2conc)
%     
%     subplot(122)
%     imagescn(pg_o2conc_filt_ds)

    
    gridded_oxy  = pg_o2conc_filt(:);
    time_oxy = pg_time(:);
    [time_oxy,idor] = sort(time_oxy);
    gridded_oxy = gridded_oxy(idor);
    [time_oxy,idu]=unique(time_oxy);
    gridded_oxy = gridded_oxy(idu);
    data.oxygen_concentration_gridded_corr(seg_idx) =interp1gap(time_oxy,gridded_oxy,timeVecSec,30);
    
    

end
data.oxygen_solubility = O2solubility(data.temperature,data.salinity_corrected);
data.oxygen_saturation = data.oxygen_concentration_gridded_corr./data.oxygen_solubility*100;
save('glider_data_oxy_processed.mat','data');


function temp = exponential_filter(var, time, tau0, tau_fact)
    % Check if tau_fact was provided
    if nargin < 3
        tau_fact = 0;
    end

    % Initialize variables
    n = length(var);
    temp = zeros(size(var));
    dt = 1; % Since time difference is always 1 second

    % Adjust tau for each point
    tau_dep = tau0 + tau_fact; 
    ef = exp(-dt ./ tau_dep);

    % Handle cases where ef is 1 or very close to 1
    ef(ef >= 0.9999) = 0.9999;

    temp(1) = var(1); % Initial condition, first point same as in original data

    for i = 2:n
        if isnan(var(i))
            temp(i) = NaN; % Retain NaNs in the output
        else
            if isnan(temp(i-1)) % If the previous value is NaN, use the current input value
                temp(i) = var(i);
            else % Regular calculation if no NaNs involved
                temp(i) = var(i) * (1 - ef) + temp(i-1) * ef;
            end
        end
    end
end

function temp = reverse_exponential_filter(var, time, tau0, tau_fact)
    % Check if tau_fact was provided
    if nargin < 3
        tau_fact = 0;
    end

    % Initialize variables
    n = length(var);
    temp = zeros(size(var));
    dt = 1; % Since time difference is always 1 second

    % Adjust tau for each point
    tau_dep = tau0 + tau_fact; 
    ef = exp(-dt ./ tau_dep);

    % Handle cases where ef is 1 or very close to 1
    ef(ef >= 0.9999) = 0.9999;

    temp(1) = var(1); % Starting condition, first point same as in original data

    for i = 2:n
        if isnan(var(i))
            temp(i) = NaN; % Retain NaNs in the output
        else
            if isnan(var(i-1)) || isnan(temp(i-1)) % If there's a NaN in previous input or output
                temp(i) = var(i); % Use the current input value
            else % Regular calculation if no NaNs involved
                temp(i) = (var(i) - var(i-1) * ef) / (1 - ef);
            end
        end
    end

end



