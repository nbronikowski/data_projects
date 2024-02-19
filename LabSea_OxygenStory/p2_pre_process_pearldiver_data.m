clear; close all; clc;

path_name = './mat_files/'
var_name  = 'pearldiver_data'
load(fullfile(path_name,[var_name,'.mat']))

%% basic QC
addpath('./processing_tools/');
% compute datetime stamps
pearldiver.time = datetime(pearldiver.time,'ConvertFrom','posixtime');

% do some basic QC'ing to get rid of weird data and too many nan's.
idbad = (pearldiver.temperature==0)|  (pearldiver.temperature==1) | (pearldiver.conductivity<1);
pearldiver.temperature(idbad)=NaN;
pearldiver.conductivity(idbad)=NaN;

% remove large conductivity spikes
idnan = find(~isnan(pearldiver.conductivity));
idSpike = abs(sgolayfilt(pearldiver.conductivity(idnan),3,11)-pearldiver.conductivity(idnan))>0.04;

figure();
hold on
plot(pearldiver.time(idnan),pearldiver.conductivity(idnan));
plot(pearldiver.time(idnan(idSpike)),pearldiver.conductivity(idnan(idSpike)),'*m')

pearldiver.conductivity(idnan(idSpike))=NaN;
pearldiver.temperature(idnan(idSpike))=NaN;
pearldiver.pressure(idnan(idSpike))=NaN;

% oxygen calphase has a lot of spikes
pearldiver.oxygen_calphase_raw = pearldiver.oxygen_calphase;
pearldiver.oxygen_calphase(pearldiver.oxygen_calphase<32)=NaN;
idnan = ~isnan(pearldiver.temperature);
pearldiver.oxygen_calphase = interp1(pearldiver.time(idnan),pearldiver.oxygen_calphase(idnan),pearldiver.time);
pearldiver.oxygen_calphase(~idnan)=NaN;

plot(pearldiver.time,pearldiver.oxygen_calphase_raw,'.b'); hold on
plot(pearldiver.time,pearldiver.oxygen_calphase,'.r'); 

yyaxis right
plot(pearldiver.time,pearldiver.temperature,'.'); 




% further step remove spikes from surface periods and science aux power
% fluctuations - SBMBII
id = isnan(pearldiver.oxygen_calphase);
temp = pearldiver.oxygen_calphase(~id);
windowSize = 15; floor(length(temp(~isnan(temp)))/1000); %
rollingStd = movstd(temp, windowSize);
rollingMean = movmean(temp, windowSize);
idSpike = abs(temp - rollingMean) > 2 * rollingStd; % k is a threshold, like 2 or 3
temp(idSpike) = NaN;
pearldiver.oxygen_calphase_cleaned = nan*pearldiver.oxygen_calphase;
pearldiver.oxygen_calphase_cleaned(~id)=temp;
clear temp

%% correct salinity - maybe best to break into profile/profile? 
idnan = find(~(isnan(pearldiver.temperature) | isnan(pearldiver.conductivity)));
conductivity=pearldiver.conductivity(idnan);
temperature=pearldiver.temperature(idnan);
pressure=pearldiver.pressure(idnan);

% Correct salinity lag
minAdv = -2.5; maxAdv = 2.5;
[bestAdv,bestSal,advances,goodness,binIndex,salinities] = find_cond_advance(...
                            temperature, conductivity, pressure,...
                            minAdv,maxAdv,'showprogress','yes');

% Apply correction to data
pearldiver.salinity_cor = NaN*pearldiver.salinity;
pearldiver.salinity_cor(idnan) = bestSal;
% figure()
% plot(pearldiver.salinity,'.'); hold on
% plot(pearldiver.salinity_cor,'.')

%% compute raw oxygen levels
% pearldiver had Oxygen Optode 4831 SN 124 calibrated last in 2018
SVUFoilCoef = [2.658710E-03,1.305085E-04,2.536294E-06,6.186505E-02,...
    5.658259E-05,-1.514249E-02,1.192067E-03];
modeltype ='uchidaAADI';
pearldiver.raw_oxygen_concentration = optcalcO2(pearldiver.temperature,...
                    pearldiver.oxygen_calphase_cleaned,SVUFoilCoef,...
                    modeltype,pearldiver.salinity_cor, 1013.25,...
                    pearldiver.pressure); 

%% Gridding
pg = 0:1:1029;%ceil(max(pearldiver.pressure,[],'omitnan')); 
timeDateNum = datenum(pearldiver.time);
[~,~,time,xu] = pgrid_columns(pearldiver.profile_index,pearldiver.pressure,timeDateNum,pg);
[~,~,oxygen_raw] = pgrid_columns(pearldiver.profile_index,pearldiver.pressure,pearldiver.raw_oxygen_concentration,pg);
[~,~,salinity] = pgrid_columns(pearldiver.profile_index,pearldiver.pressure,pearldiver.salinity_cor,pg);
[~,~,temperature] = pgrid_columns(pearldiver.profile_index,pearldiver.pressure,pearldiver.temperature,pg);
[~,~,conductivity] = pgrid_columns(pearldiver.profile_index,pearldiver.pressure,pearldiver.conductivity,pg);
[~,~,pressure] = pgrid_columns(pearldiver.profile_index,pearldiver.pressure,pearldiver.pressure,pg);

[pearldiver.gridded.oxygen_raw,column_idx]=deleteAlmostEmptyColumns(oxygen_raw,pg);
% [pearldiver.gridded.temperature,~]=deleteAlmostEmptyColumns(temperature,pg);
% [pearldiver.gridded.salinity,~]=deleteAlmostEmptyColumns(salinity,pg);
% [pearldiver.gridded.conductivity,~]=deleteAlmostEmptyColumns(conductivity,pg);
pearldiver.gridded.time = time(:,column_idx);
pearldiver.gridded.profile_index = xu(column_idx);
pearldiver.gridded.pressure_grid = pg;
% pearldiver.gridded.pressure = deleteAlmostEmptyColumns(pressure,pg);
pearldiver.gridded.temperature = temperature(:,column_idx);
pearldiver.gridded.salinity = salinity(:,column_idx);
pearldiver.gridded.conductivity = conductivity(:,column_idx);
pearldiver.gridded.pressure = pressure(:,column_idx);

%% Add extra vars
pearldiver.dateNum = datenum(pearldiver.time);
pearldiver.gridded.timeg = mean(pearldiver.gridded.time,1,'omitnan');
pearldiver.gridded.lon = interp1(pearldiver.dateNum,pearldiver.longitude,pearldiver.gridded.timeg);
pearldiver.gridded.lat = interp1(pearldiver.dateNum,pearldiver.latitude,pearldiver.gridded.timeg);

save(fullfile(path_name,[var_name,'_clean.mat']),'pearldiver')
