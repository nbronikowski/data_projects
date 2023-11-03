clear; close all; clc;

path_name = './mat_files/'
var_name  = 'sunfish_data.mat'
load(fullfile(path_name,var_name))

%% Basic Cleaning (not really scientific)

% compute datetime stamps
sunfish.time = datetime(sunfish.time,'ConvertFrom','posixtime');

% do some basic QC'ing to get rid of weird data and too many nan's.
idbad = (sunfish.temperature==1) | (sunfish.conductivity<1);
sunfish.temperature(idbad)=NaN;
sunfish.conductivity(idbad)=NaN;

% process oxygen
% sunfish had optode 673
SVUFoilCoef = [2.78565E-03,1.17409E-04,2.39332E-06,2.32859E+02,-3.80398E-01,-4.79825E+01,4.62200E+00];
modeltype ='uchidaAADI';
sunfish.raw_oxygen_concentration = optcalcO2(sunfish.temperature,sunfish.oxygen_calphase,SVUFoilCoef, modeltype,sunfish.salinity, 1013.25,sunfish.pressure); 

save(fullfile(path_name,[var_name,'_oxy.mat']),'sunfish')


