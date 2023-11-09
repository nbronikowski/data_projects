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


%% Gridding
pg = 0:1:ceil(max(sunfish.pressure,[],'omitnan')); 
timeDateNum = datenum(sunfish.time);
[~,~,time,xu] = pgrid_columns(sunfish.profile_index,sunfish.pressure,timeDateNum,pg);
[~,~,oxygen_raw] = pgrid_columns(sunfish.profile_index,sunfish.pressure,sunfish.raw_oxygen_concentration,pg);
[~,~,salinity] = pgrid_columns(sunfish.profile_index,sunfish.pressure,sunfish.salinity_cor,pg);
[~,~,temperature] = pgrid_columns(sunfish.profile_index,sunfish.pressure,sunfish.temperature,pg);
[~,~,conductivity] = pgrid_columns(sunfish.profile_index,sunfish.pressure,sunfish.conductivity,pg);
[~,~,pressure] = pgrid_columns(sunfish.profile_index,sunfish.pressure,sunfish.pressure,pg);

[sunfish.gridded.temperature,~]=deleteAlmostEmptyColumns(temperature,pg);
[sunfish.gridded.salinity,~]=deleteAlmostEmptyColumns(salinity,pg);
[sunfish.gridded.conductivity,column_idx]=deleteAlmostEmptyColumns(conductivity,pg);
sunfish.gridded.time = time(:,column_idx);
sunfish.gridded.profile_index = xu(column_idx);
sunfish.gridded.pressure_grid = pg;
sunfish.gridded.pressure = deleteAlmostEmptyColumns(pressure,pg);
sunfish.gridded.oxygen_raw = oxygen_raw(:,column_idx);


%% Save result
save(fullfile(path_name,[var_name,'_clean.mat']),'sunfish')


