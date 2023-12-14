clear; close all; clc;

% Update Nov 9
% - missing ecopuck processing for now

path_name = './mat_files/'
var_name  = 'sunfish_data'
load(fullfile([path_name,var_name,'.mat']))

%% basic QC

% compute datetime stamps
sunfish.time = datetime(sunfish.time,'ConvertFrom','posixtime');

% do some basic QC'ing to get rid of weird data and too many nan's.
idbad = (sunfish.temperature==1) | (sunfish.conductivity<1);
sunfish.temperature(idbad)=NaN;
sunfish.conductivity(idbad)=NaN;

% remove large conductivity spikes
idnan = find(~isnan(sunfish.conductivity));
idSpike = abs(sgolayfilt(sunfish.conductivity(idnan),3,11)-sunfish.conductivity(idnan))>0.04;

% figure();
% hold on
% plot(sunfish.time(idnan),sunfish.conductivity(idnan));
% plot(sunfish.time(idnan(idSpike)),sunfish.conductivity(idnan(idSpike)),'*m')

sunfish.conductivity(idnan(idSpike))=NaN;
sunfish.temperature(idnan(idSpike))=NaN;
sunfish.pressure(idnan(idSpike))=NaN;

%% correct salinity - maybe best to break into profile/profile? 
idnan = find(~(isnan(sunfish.temperature) | isnan(sunfish.conductivity)));
conductivity=sunfish.conductivity(idnan);
temperature=sunfish.temperature(idnan);
pressure=sunfish.pressure(idnan);

% Correct salinity lag
minAdv = -2.5; maxAdv = 2.5;
[bestAdv,bestSal,advances,goodness,binIndex,salinities] = find_cond_advance(...
                            temperature, conductivity, pressure,...
                            minAdv,maxAdv,'showprogress','yes');
% Apply correction to data
sunfish.salinity_cor = NaN*sunfish.salinity;
sunfish.salinity_cor(idnan) = bestSal;
figure()
plot(sunfish.salinity,'.'); hold on
plot(sunfish.salinity_cor,'.')


%% Compute raw oxygen
% sunfish had oxygen optode SN 673 with factory calibration from 2017
SVUFoilCoef = [2.78565E-03,1.17409E-04,2.39332E-06,2.32859E+02,-3.80398E-01,-4.79825E+01,4.62200E+00];
modeltype ='uchidaAADI';
sunfish.raw_oxygen_concentration = optcalcO2(sunfish.temperature,...
                sunfish.oxygen_calphase,SVUFoilCoef, modeltype,...
                sunfish.salinity_cor, 1013.25,sunfish.pressure); 


%% Ecopuck Processing
addpath('./processing_tools/');

% Dark Counts Calc. in-situ
sunfish.dateTime = sunfish.time;
sunfish.dateNum = datenum(sunfish.time);
data_range = sunfish.pressure>350 & sunfish.pressure<900 & ...
    sunfish.dateNum>datenum(2022,01,20) & sunfish.dateNum<datenum(2022,04,01);

% Backscatter 700nm 
eco_bb_dark_counts = 47; % factory value
eco_bb_scale_factor = 1.913E-06; % m-1 sr-1 / counts
% eco_bb_dark_counts = 110; %min_75th_percentile(sunfish.bb700_sig(sunfisha_range))

lambda = 700; % 700 nm backscatter 
theta = 117; % angle of instrument
VBSC  = eco_bb_scale_factor*(sunfish.bb700_sig-eco_bb_dark_counts);  

% 1. Method Zhang et al., 2009 Volume Scattering of Seawater
beta_sw=betasw(lambda, sunfish.temperature, theta, sunfish.salinity_cor);

% Particle backscatter coefficient, BBP 
Xp = 1.1; %recommended by WETLabs 2013
sunfish.bbp700 = 2*pi*(VBSC - beta_sw) * Xp;

% Chlorophyll
eco_chl_dark_counts = 45; % factory value
eco_chl_scale_factor = 0.0073; % micro-g L-1
eco_chl_dark_counts = min_75th_percentile(sunfish.chlor_sig(data_range))
sunfish.chlor  = eco_chl_scale_factor*(sunfish.chlor_sig-eco_chl_dark_counts);  

% CDOM
eco_cdom_dark_counts = 50;
eco_cdom_scale_factor = 0.0904;
eco_cdom_dark_counts = min_75th_percentile(sunfish.cdom_sig(data_range))
sunfish.cdom  = eco_cdom_scale_factor*(sunfish.cdom_sig-eco_cdom_dark_counts);  

% remove weird data
sunfish.chlor(sunfish.chlor<0)=NaN;
sunfish.bbp700(sunfish.bbp700<0)=NaN;
sunfish.cdom(sunfish.cdom<0)=NaN;

%% Gridding
pg = 0:1:ceil(max(sunfish.pressure,[],'omitnan')); 
timeDateNum = datenum(sunfish.time);
[~,~,time,xu] = pgrid_columns(sunfish.profile_index,sunfish.pressure,timeDateNum,pg);
[~,~,oxygen_raw] = pgrid_columns(sunfish.profile_index,sunfish.pressure,sunfish.raw_oxygen_concentration,pg);
[~,~,salinity] = pgrid_columns(sunfish.profile_index,sunfish.pressure,sunfish.salinity_cor,pg);
[~,~,temperature] = pgrid_columns(sunfish.profile_index,sunfish.pressure,sunfish.temperature,pg);
[~,~,conductivity] = pgrid_columns(sunfish.profile_index,sunfish.pressure,sunfish.conductivity,pg);
[~,~,pressure] = pgrid_columns(sunfish.profile_index,sunfish.pressure,sunfish.pressure,pg);

[~,~,chlor] = pgrid_columns(sunfish.profile_index,sunfish.pressure,sunfish.chlor,pg);
[~,~,cdom] = pgrid_columns(sunfish.profile_index,sunfish.pressure,sunfish.cdom,pg);
[~,~,bbp700] = pgrid_columns(sunfish.profile_index,sunfish.pressure,sunfish.bbp700,pg);

[sunfish.gridded.temperature,~]=deleteAlmostEmptyColumns(temperature,pg);
[sunfish.gridded.salinity,~]=deleteAlmostEmptyColumns(salinity,pg);
[sunfish.gridded.conductivity,column_idx]=deleteAlmostEmptyColumns(conductivity,pg);
sunfish.gridded.time = time(:,column_idx);
sunfish.gridded.profile_index = xu(column_idx);
sunfish.gridded.pressure_grid = pg;
sunfish.gridded.pressure = deleteAlmostEmptyColumns(pressure,pg);
sunfish.gridded.oxygen_raw = oxygen_raw(:,column_idx);
sunfish.gridded.chlor = chlor(:,column_idx);
sunfish.gridded.cdom  = cdom(:,column_idx);
sunfish.gridded.bbp700 = bbp700(:,column_idx);

%% Add extra vars
sunfish.gridded.timeg = mean(sunfish.gridded.time,1,'omitnan');
sunfish.gridded.lon = interp1(sunfish.dateNum,sunfish.longitude,sunfish.gridded.timeg);
sunfish.gridded.lat = interp1(sunfish.dateNum,sunfish.latitude,sunfish.gridded.timeg);

%% Save result
save(fullfile(path_name,[var_name,'_clean.mat']),'sunfish','-v7.3')
