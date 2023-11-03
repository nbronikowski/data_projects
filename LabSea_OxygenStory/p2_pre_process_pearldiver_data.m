clear; close all; clc;

path_name = './mat_files/'
var_name  = 'pearldiver_data'
load(fullfile(path_name,[var_name,'.mat']))

%% Basic Cleaning (not really scientific)

% compute datetime stamps
pearldiver.time = datetime(pearldiver.time,'ConvertFrom','posixtime');

% do some basic QC'ing to get rid of weird data and too many nan's.
idbad = (pearldiver.temperature<1) | (pearldiver.conductivity<1);
pearldiver.temperature(idbad)=NaN;
pearldiver.conductivity(idbad)=NaN;

% oxygen calphase has a lot of spikes
idbad = (pearldiver.oxygen_calphase<33) | isnan(pearldiver.temperature);
pearldiver.oxygen_calphase(idbad)=NaN;

% further step remove spikes from surface periods and science aux power
% fluctuations - SBMBII
id = isnan(pearldiver.oxygen_calphase);
temp = pearldiver.oxygen_calphase(~id);
windowSize = 15; floor(length(temp(~isnan(temp)))/1000); %
rollingStd = movstd(temp, windowSize);
rollingMean = movmean(temp, windowSize);
idSpikes = abs(temp - rollingMean) > 2 * rollingStd; % k is a threshold, like 2 or 3
temp(idSpikes) = NaN;
pearldiver.oxygen_calphase_cleaned = nan*pearldiver.oxygen_calphase;
pearldiver.oxygen_calphase_cleaned(~id)=temp;

% compute raw oxygen levels
% pearldiver had optode 124
SVUFoilCoef = [2.658710E-03,1.305085E-04,2.536294E-06,6.186505E-02,...
    5.658259E-05,-1.514249E-02,1.192067E-03];
modeltype ='uchidaAADI';
pearldiver.raw_oxygen_concentration = optcalcO2(pearldiver.temperature,pearldiver.oxygen_calphase_cleaned,SVUFoilCoef, modeltype,pearldiver.salinity, 1013.25,pearldiver.pressure); 

save(fullfile(path_name,[var_name,'_oxy.mat']),'pearldiver')

%% Oxygen Correction - trying the GEOMAR way



%% Gridding