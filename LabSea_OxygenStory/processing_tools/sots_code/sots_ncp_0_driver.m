clear all; clc; close all; clear;
% Southern Ocean Time Series Net Community Production driver script
% B Weeding and T W Trull 2020

% This script reads input data and calls other scripts to calculate NCP from times series observations.


% The input data must include time series of the following variables, or default choices: 
% 1. mooring_data.time = the times of the observations in format:
% 2. mooring_data.gastension_Pa = GTD pressure in Pa
% 3. mooring_data.dox2_umolkg = Dissolved Oxygen concentration in umol/kg
% 4. mooring_data.mld_m = mixed layer depth in meters
% 5. mooring_data.temp_C = mixed layer temperature in C
% 6. mooring_data.psal_PSU = mixed layer salinity in psu
% 7. mooring_data.windspeed_ms = wind speed at 10m height in m s-1
% 8. mooring_data.atmosphericpress_Pa = atmospheric pressure at sea level in Pa 
% (a default value can be set if not available) 
% 9. mooring_data.dox2_submld_umolkg = subsurface oxygen concentration
% (a default value can be set if not available)


% The script is designed to read input data from SOTS IMOS netcdf format files, and place it in a Matlab structure called mooring_data. 
% To use with other data, create the mooring_data structure using the variable names as listed above.

% The script requires that the input data and the following scripts must all be in the Matlab path:
% 1. sots_ncp_1_addpath				% modify to add and save the path specific to your computer
% 2. sots_ncp_2_constants		% loads scientific constants
% 3. sots_ncp_3_extractor		% loads data from SOTS netCDF files into structure mooring_data
% 4. sots_ncp_4_select_acceptable_qc	% selects data using flags in SOTS file, looking for the longest continuous period of acceptably flagged data
% 5. sots_ncp_5_mooring_plot		% makes plots to check that data is present and correct
% 6. sots_ncp_6_dox2_converter    % convertss O2 and solubility from umol/kg to mol/m^3, calculates saturation un umol/kg
% 7. sots_ncp_7_dox2_sat_plot     % produces plots of O2 saturation with temperature and gas tension
% 8. sots_ncp_8_gas_records       % calculates Henry's law constants, water vapour and N2 time series
% 9. sots_ncp_9_gas_transfer_velocity_X   % calculates gas transfer velocity according to user choice
% 10. sots_ncp_10_N2_exchange      % calculates air-sea exchange and bubble processes to match N2 time series
% 11. sots_ncp_11_O2_exchange_X     % calculates O2 exchange processes to produce physical O2 timeseries
% 12. sots_ncp_12_atmospheric_exchange_visualisation   %creates plots to compare N2 and physical O2 exchange processes  
% 13. sots_ncp_13_net_com_prod_calc    % calculates net community production
% 14. sots_ncp_14_net_com_prod_visualisation   % creates plots to compare NCP time series with relevant variables
% 15. sots_ncp_15_workspace_cleanup     % deletes unwanted variables and assigns the variables described in 'Outputs'

% all of these are required except: 5, 7, 12, 14, 15

% The script allows choices for different algorithms and parameters, as follows.

%% Choice 0. Select the deployment
% Currently available options are 'pulse7','pulse9','sofs75'

% deployment_choice = 'pulse9';

%% Choice 1. Set the air-sea gas exchange algorithm: 
airsea_algorithm = 1;  

%1 is the default, the full set of options are:
% 1. Wanninkhof 2014
% 2. Wanninkhof 2002
% 3. Ho 2006
% 4. Sweeney 2007
% 5. Edson 2011

%% Choice 2. Set the parameter, bubble_beta, which is the ratio of fully to partially dissolving bubbles 
bubble_beta = 10;

% 1 is the default, values of 0.1 and 10 provide useful bounds to explore 
% the sensitivity of NCP to this unknown aspect of the fate of injected 
% bubbles. See Emerson et al. 2008 (eqn 7) for further information. 

%% Choice 3. Set the vertical eddy diffusivity and gradient thickness

eddy_diff_coeff = 0.33E-4;

% The default eddy diffusivity coefficient is 0.33x10^-4 in m^2 s^-1 from
% Law et al. 2003 .

eddy_gradient_thickness_m = 50;

% The default value is 50m from Weeding and Trull 2014, based on CTD and Argo O2 profiles at SOTS
% if a subsurface O2 record is used, this depth should be evaluated to match the sensor spacing and associated O2 gradients

%% Choice 4. Set the time span over which the MLD estimates are smoothed prior to use 

mld_smooth_span = 72;

% The default value is 24 hours, to emphasize the longer time scales on 
% which the mixed layer and sub-mixed layer water masses are likely to have 
% different oxygen contents

%% Choice 5. Set the default atmospheric pressure if the variable is not available
%If there is an atmospheric pressure record, it is used by default.
%To set a constant atmospheric pressure, set:
%atmospheric_pressure_manual_overide=1
%and your choice of atmospheric pressure in atmospheres, e.g. for 1 atm=101325Pa, set:
%atmospheric_pressure_choice = 1
%[the default value of atmospheric_pressure_choice is 0]

atmospheric_pressure_choice = 1;

% This sets a steady atmospheric pressure value in atm(101325Pa), if a timeseries is
% not available. The presence of a timeseries
% (mooring_data.atmosphericpress_Pa) will override this value.

atmospheric_pressure_manual_override = 0; % default 0

% Setting 'atmospheric_pressure_manual_override' to 1 will override the available 
% atmospheric pressure vector from the mooring, and use the user specified constant
% as set by atmospheric_pressure_choice (in atm).

%% Choice 6. Set the default subsurface oxygen concentration if the variable is not available
%If there is a subsurface oxygen concentration, it is used by default.
%To set a constant subsurface oxygen concentration, set:
%sub_mld_dox2_manual_override =1
%and your choice of subsurface oxygen concentration in moles m-3, e.g. for 264 uM oxygen, set:
%sub_mld_dox2_choice = 0.264;
%[the default value of sub_mld_dox2_choice is 0]
sub_mld_dox2_choice = 0;

sub_mld_dox2_choice = 0.28;

% This sets a steady oxygen value below the mixed layer in mol/m^3, if a
% timeseries is not available. The presence of 'mooring_data.sub_mld_dox2_umolkg'
% will override any value given here. 0.264 (as used in Weeding and Trull, 2014) 
% at a density of 1029.15kg/m^3 (typical 500m density) corresponds to ~256.5umol/kg
% An alternate useful choice is the first value in a timeseries, especially
% when deep mixing has occured eg. 0.28 for Pulse9.
 
sub_mld_dox2_manual_override = 0; % default 0

% Setting 'sub_mld_dox2_manual_override' to 1 will override the available 
% sub MLD O2 vector from the mooring, and use the user specified constant
% as set by sub_mld_dox2_choice (in mol/m^3).


%% Choice 7. Decide whether to include eddy diffusion andentrainment

exchange_choice = 1;

% 1: only air-sea exchange and bubbles
% 2: option 1 plus eddy diffusion 
% 3: option 2 plus MLD based entrainment 
% 4: option 2 plus and heat budget based entrainment

%% Choice 8. Selecting qc limits
% Set this parameter to the lowest qc value you're not willing to accept.
% Default value is 3, so that only data of qc values 1 and 2 are used by 
% sots_ncp_select_acceptable_qc.
qc_limit = 3;


%% Outputs
% The scripts create outputs in the Matlab workspace.  The important ones to save from a run are all the choices (as listed above), and the calculated results:
% 1. NCP					% hourly timeseries of net community production in moles of O2 per cubic metre 
% 2. NCP_qc                 % qc values for NCP, calculated as the maximum qc value of the mooring_data variables
% 3. Bubble_injected_oxygen	% hourly timeseries of bubble injected oxygen in moles per cubic metre
% 4. Entrained_oxygen		% hourly timeseries of entrained oxygen (from below MLD) in moles per cubic metre
% 5. Eddy_diff_oxygen       % hourly timeseries of eddy diffused oxygen (from below MLD) in moles per cubic metre
% 6. Gas_exchange_oxygen    % hourly timeseries of oxygen exchanged into water from atmosphere in moles per cubic metre

% Cumulatively summed timeseries of NCP per square metre are available in the mooring_data structure
% ie. mooring_data.ncp_O2_umm2_cumsum and mooring_data.ncp_C_mgm2_cumsum
%% Run the scripts

% sots_ncp_1_addpath
 
sots_ncp_2_constants; 

sots_ncp_3_extractor; 

% sots_ncp_4_select_acceptable_qc

sots_ncp_5_mooring_plot

sots_ncp_6_O2_converter 

sots_ncp_7_O2_sat_plot

sots_ncp_8_gas_records 

switch airsea_algorithm
    case 1 % Wanninkhof 2014
        sots_ncp_9_gas_transfer_velocity_W14
        
    case 2 % Wanninkhof 2002
        sots_ncp_9_gas_transfer_velocity_W02
        
    case 3 % Ho 2006
        sots_ncp_9_gas_transfer_velocity_W02
        
    case 4 % Sweeney 2007
        sots_ncp_9_gas_transfer_velocity_S07
        
    case 5 % Edson 2011
        sots_ncp_9_gas_transfer_velocity_E11
end

sots_ncp_10_N2_exchange 

switch exchange_choice
    
    case 1 % No entrainment or eddy diffusion
        sots_ncp_11_O2_exchange
        
    case 2 % Eddy diffusion no entrainment
        sots_ncp_11_O2_exchange_eddy
        
    case 3 % Eddy diffusion and entrainment from MLD record
        sots_ncp_11_O2_exchange_eddy_and_entrainment
        
    case 4 % Entrainment from heat budget MLD

end

sots_ncp_12_atmospheric_exchange_visualisation

sots_ncp_13_net_com_prod_calc 

sots_ncp_14_net_com_prod_visualisation

% sots_ncp_15_workspace_cleanup

