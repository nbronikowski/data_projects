clear; close all; clc;


load('glider_data.mat');
gdat.time = gdat.timeDateNum;


idnan = ~isnan(gdat.pres);
[prof_idx,prof_dir] = findProfiles(gdat.time(idnan),gdat.pres(idnan),'STALL',20);
gdat.prof_idx(~idnan) = NaN;
gdat.prof_idx(idnan)=prof_idx;
gdat.prof_dir(idnan)=prof_dir;
gdat.prof_dir(~idnan)=NaN;
%% Segement: we extract every up /down cast
timeVec = gdat.time;
uprofs = unique(gdat.prof_idx(rem(gdat.prof_idx,1)==0));
k = 1;
for i = 1:length(uprofs)-1
    idx1=find(gdat.prof_idx == uprofs(i)   & gdat.prof_dir == 1   ); 
    idx2=find(gdat.prof_idx == uprofs(i+1) & gdat.prof_dir == -1  ); 
    
    if length(idx1)>10 && length(idx2)>10
        idx = [idx1;idx2];
        deltaTime = (max(timeVec(idx),[],'omitnan')...
            -min(timeVec(idx),[],'omitnan'))*24*60; % mins
        if deltaTime>5
            gdat.SegStartTime(k) = min(timeVec(idx),[],'omitnan');
            gdat.SegEndTime(k) = max(timeVec(idx),[],'omitnan');
            k = k + 1;
%             figure(1)
%             plot(timeVec(idx),-gdat.pres(idx),'.'); hold on
%             cla
        end
    else
    end
    
end

%% Segment Finding for finding true segments (between surfacings)
% Surface Times x_dr_state > 1, 2 = just surfaced, 4=about to dive

% drState = gdat.x_dr_state;
% tempIDX=find(gdat.x_dr_state>1); % surface positions
% i_dive = tempIDX(diff(tempIDX)>1);
% i_sfce = tempIDX(([0;diff(tempIDX)])>1);
% 
% diveTimes = timeVec(i_dive);
% surfTimes = timeVec(i_sfce);
% 
% k = 1;
% for i = 1:length(diveTimes)
%     if (surfTimes(i)-diveTimes(i))*24*60>20 % at least 30 min segment
%         gdat.SegStartTime(k)=diveTimes(i);
%         gdat.SegEndTime(k)=surfTimes(i);
%         k = k+1;
%         % SegTimeIdx=find(timeVec>=diveTimes(i) & timeVec<=surfTimes(i)); 
%         % plot(timeVec(SegTimeIdx),-gdat.m_depth(SegTimeIdx),'.'); hold on
%     end
% end

%% Grab Nortek data using parse_nortek_adcp script
dat = parse_nortek_adcp('./adcp_data/Data.ad2cp');
adat.time = dat.nortek_avg.time;

% combine data and make sure its over the time period that matches the
% glider data. ADCP grabs time from glider so it should be mostly in sync
idx = find(adat.time>=nanmin(gdat.time) & adat.time<=nanmax(gdat.time));

adat.time = adat.time(idx);
adat.pres = dat.nortek_avg.pressure(idx); % pressure
adat.roll = dat.nortek_avg.roll(idx);
adat.heading = dat.nortek_avg.heading(idx);
adat.pitch = dat.nortek_avg.pitch(idx);
adat.ADCP_sos  = dat.nortek_avg.sos(idx); % assuming this is the speed of sound used by Nortek instead of 1500 m/s??
adat.gliderSegStartTime = gdat.SegStartTime';
adat.gliderSegEndTime = gdat.SegEndTime';

adat.Vm = dat.nortek_avg.vel(:,:,idx);
adat.Am = dat.nortek_avg.amp(:,:,idx);
adat.Cr = dat.nortek_avg.cor(:,:,idx);
adat.Magx = dat.nortek_avg.raw_mag(idx,:);
adat.bt_range = dat.bt_range(:,idx);
adat.bt_vel = dat.bt_vel(:,idx);
adat.bt_time = dat.bt_time(idx); 

adat.Acc = dat.nortek_avg.raw_acc(idx,:);
adat.water_temperature = dat.nortek_avg.temp(idx);
adat.cell_size = dat.nortek_avg.cell_size(idx);
adat.cor_nom = dat.nortek_avg.cor_nom(idx);
adat.bat_vol = dat.nortek_avg.bat_voltage(idx);
adat.blanking = dat.nortek_avg.blanking(idx);
adat.pres_temp = dat.nortek_avg.pres_temp(idx);
adat.ncells = dat.nortek_avg.ncells;
adat.nbeams = dat.nortek_avg.nbeams;
adat.serial = dat.nortek_avg.serial;

% FROM HDR File:
adat.beam2xyz = ...
        [0.6782	    0	   -0.6782	    0      ;...
        0	   -1.1831	    0	        1.1831 ;...
        0.74	    0	    0.74	    0      ;...
        0	    0.5518	    0	        0.5518];

adat.range_cells = (1:adat.ncells)+mean(adat.blanking);


%% ADD GLIDER DATA VARIABLES
% glider RBR CTD pressure is more accurate
idnan = ~isnan(gdat.pres);
adat.glider_pres = interp1gap(gdat.time(idnan),...
    gdat.pres(idnan),adat.time,300/86400);

idnan = ~isnan(gdat.depth);
adat.glider_depth = interp1gap(gdat.time(idnan),...
    gdat.depth(idnan),adat.time,300/86400);

idnan = ~isnan(gdat.c_stp);
adat.glider_CTD_sos = interp1gap(gdat.time(idnan),...
    gdat.c_stp(idnan),adat.time,300/86400);

idnan = ~isnan(gdat.temp);
adat.glider_CTD_temp = interp1gap(gdat.time(idnan),...
    gdat.temp(idnan),adat.time,30/86400);

idnan = ~isnan(gdat.salt);
adat.glider_CTD_salt = interp1gap(gdat.time(idnan),...
    gdat.salt(idnan),adat.time,30/86400);

% roll, heading, pitch
idnan = ~isnan(gdat.m_pitch);
adat.glider_pitch = interp1(gdat.time(idnan),...
    gdat.m_pitch(idnan),adat.time);

idnan = ~isnan(gdat.m_heading);
adat.glider_heading = interp1(gdat.time(idnan),...
    gdat.m_heading(idnan),adat.time);

idnan = ~isnan(gdat.m_roll);
adat.glider_roll = interp1(gdat.time(idnan),...
    gdat.m_roll(idnan),adat.time);

% water depth
idnan = ~isnan(gdat.water_depth);
adat.glider_water_depth = interp1(gdat.time(idnan),...
    gdat.water_depth(idnan),adat.time);

% long & lat & u, v
idnan = ~isnan(gdat.lon);
adat.glider_lon = interp1(gdat.time(idnan),...
    gdat.lon(idnan),adat.time);

idnan = ~isnan(gdat.lat);
adat.glider_lat = interp1(gdat.time(idnan),...
    gdat.lat(idnan),adat.time);

% Interpolated glider currents to match segments
idnan = ~isnan(gdat.u_dac);
adat.glider_u_dac = interp1(gdat.time(idnan),...
    gdat.u_dac(idnan),adat.time);

idnan = ~isnan(gdat.v_dac);
adat.glider_v_dac = interp1(gdat.time(idnan),...
    gdat.v_dac(idnan),adat.time);

% Also keep a copy of not interpolated u,v DACs
adat.glider_u_dac_o = adat.glider_u_dac*NaN;
adat.glider_v_dac_o = adat.glider_v_dac*NaN;
uv_nnan_times = gdat.time(idnan);

% For each time in the new array, check if it matches a NaN time from the original data. If it does, set the new data to NaN
for i = 1:length(uv_nnan_times)
    [~,idx]=min(abs(adat.time-uv_nnan_times(i)),[],'omitnan');
    adat.glider_u_dac_o(idx) = adat.glider_u_dac(idx);
    adat.glider_v_dac_o(idx) = adat.glider_v_dac(idx);
end

save('ADCP+glider_data.mat','adat')
