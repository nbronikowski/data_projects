clear; clc; close all;
load unit_334_delayed_trajectory.nc.mat

%% Remove Bad data and Interpolate Data

gdat.cond(gdat.cond<1)=NaN;
gdat.salt(gdat.salt<10)=NaN;
gdat.temp(gdat.temp==0)=NaN;

idnan = ~isnan(gdat.temp) | ~isnan(gdat.salt) | ~isnan(gdat.pres) ;
gdat.temp=interp1gap(gdat.timeDateNum(idnan),gdat.temp(idnan),gdat.timeDateNum,30/86400);
gdat.salt=interp1gap(gdat.timeDateNum(idnan),gdat.salt(idnan),gdat.timeDateNum,30/86400);
gdat.pres=interp1gap(gdat.timeDateNum(idnan),gdat.pres(idnan),gdat.timeDateNum,30/86400);

%% Salinity Correction
fn = fieldnames(gdat);
timeVec = gdat.timeDateNum;
data2=gdat;
for i = 1:length(fn)
    temp = data2.(fn{i});
    idnan = find(~isnan(temp)); % indices where there are non-NaN values
    if ~isempty(idnan) % check if there are any non-NaN values
        % Now, find where the NaNs are
        id_is_nan = find(isnan(temp)); % indices where there are NaN values
        if ~isempty(id_is_nan) % ensure there are NaNs to replace
            % Perform interpolation
            temp(id_is_nan) = interp1(timeVec(idnan), temp(idnan), ...
                timeVec(id_is_nan), 'linear','extrap');
            data2.(fn{i}) = temp;
        end
    else
        data2.(fn{i}) = NaN(size(data2.(fn{i})));
    end
end

gdat.salt_corr = NaN*gdat.salt;
idSeg = [1;find([0;diff(gdat.timeDateNum)]>0.5);length(gdat.timeDateNum)];
for i = 1:length(idSeg)-1

    % deal with salinity
    seg_idx = idSeg(i):idSeg(i+1)-1; % gives the section index
    tempNaN = data2.temp(seg_idx);
    condNaN = data2.cond(seg_idx);
    presNaN = data2.pres(seg_idx);
    idNaN = find(~isnan(condNaN));
    tempNaN = tempNaN(idNaN);
    condNaN = condNaN(idNaN);
    presNaN = presNaN(idNaN);
    minAdv = -10; maxAdv = 10;
    [bestAdv,bestSal,~,~,~,~] = find_cond_advance(tempNaN,condNaN,presNaN,...
                    minAdv,maxAdv,'showprogress','yes');    
    gdat.salt_corr(seg_idx)=bestSal(idNaN);
end

save('glider_data_CTD_corrected.mat','gdat');

