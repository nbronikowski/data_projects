

% [gdat.prof_idx,gdat.profile_direction] = findProfiles(gdat.time, ...
%         gdat.sci_water_pressure,'STALL',20);

uidx = unique(gdat.prof_idx);
uidx = uidx(mod(uidx,1)==0);
idx = (mod(gdat.prof_idx,1)==0);
gdat.prof_idx(~idx) =NaN;

for i = 1:length(uidx)
    idx = find(gdat.prof_idx == uidx(i));
    
    els = gdat.prof_idx(idx);
    if length(els(~isnan(els)))<10
        gdat.prof_idx(idx)=NaN;
        gdat.prof_idx(idx);
    end
end

% plot(gdat.prof_idx,gdat.oxygen_concentration,'.')


% gridding routine based on u_prof_idx
pg = 0:1:340; 

% pgrid in coloums

[pg_salt_nan,~,pg_salt] = pgrid_columns(gdat.prof_idx,gdat.depth,gdat.salt_corrected,pg);
[pg_dens0_nan,~,pg_dens0] = pgrid_columns(gdat.prof_idx,gdat.depth,gsw_sigma0(gdat.saltSA,gdat.ctemp),pg);
[pg_temp_nan,~,pg_temp,xu] = pgrid_columns(gdat.prof_idx,gdat.depth,gdat.temp,pg);
[pg_time_nan,~,pg_time,~] = pgrid_columns(gdat.prof_idx,gdat.depth,gdat.timeDateNum,pg);
[pg_o2conc_nan,~,pg_o2conc,~] = pgrid_columns(gdat.prof_idx,gdat.depth,gdat.oxygen_phase_recalc,pg);

timeg  = nanmean(pg_time,1);
D1 = 1:length(timeg); idn = ~isnan(timeg);
timeg  = interp1(D1(idn),timeg(idn),D1);

[pg_o2conc,idx]=deleteAlmostEmptyColumns(pg_o2conc,pg);
pg_temp = pg_temp(:,idx);
pg_salt = pg_salt(:,idx);
pg_time = pg_time(:,idx);
pg_dens0 = pg_dens0(:,idx);
timeg = timeg(idx);
[x,y] = meshgrid(timeg,pg);
imagescn(x,-y,CORR_DOXY_V2')
% % % 


% % % % [ tau ] = calculate_tau(pg_time',pg_pres',pg_o2conc','zres',0.25,'tlim',[0,150]);
% % % % [CORR_DOXY] = apply_tau(pg_time',pg_pres',pg_o2conc',medfilt1(tau,25));

% [thickness, tau_Tref] = calculate_tau_wTemp(pg_time',y',pg_o2conc',pg_temp');
load('./glider_data/tau_result.mat');
[CORR_DOXY] = apply_tau_wTemp(pg_time',y',pg_o2conc',pg_temp',thickness);
[CORR_DOXY_V2] = apply_tau_wTemp(pg_time',y',pg_o2conc',pg_temp',medfilt1(thickness,20));

%%
oxy= medfilt2(CORR_DOXY_V2',[2 2],'symmetric');
oxy = oxy(:);
time_oxy = pg_time(:);
[time_oxy,idor] = sort(time_oxy);
oxy = oxy(idor);

[time_oxy,idu]=unique(time_oxy);
oxy = oxy(idu);
gdat.oxygen_corrected = interp1gap(time_oxy,oxy,gdat.timeDateNum,100/86400);


gdat.gridded.oxy_corr=medfilt2(CORR_DOXY_V2',[2 2]);
gdat.gridded.oxy = pg_o2conc;
gdat.gridded.salt = pg_salt;
gdat.gridded.dens0 = pg_dens0;
gdat.gridded.temp = pg_temp;
gdat.gridded.time = timeg;
gdat.gridded.depthg = pg;
gdat.gridded.prof_idx = xu;
gdat.gridded.time_gridded = pg_time;

save(['glider_data_oxy.mat'],'gdat','-v7.3')

