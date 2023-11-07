clear; clc; close all;

%% Load Data
load ADCP+glider_data.mat                    

%% Correct Magnetometer readings
% Copying the code from Sam Coakley 
% https://github.com/sam12396/AD2CP_processing/blob/master/AD2CP_mag_correct.m
ptcheck=0;
adat.heading_corr=correct_AD2CP_heading(adat.time,adat.heading,adat.pitch,...
    adat.roll,adat.pres,adat.Magx(:,1),adat.Magx(:,2),adat.Magx(:,3),ptcheck);
if ptcheck>1
    save_figure(gcf,['./plots/ad2cp_heading_corr'],[7.5 5],'.png','300')
    close gcf
end
%% Correction Visbeck 2002 - local speed of sound correction.
% Using speed of sound from ADCP that NORTEK gave me instead of 1500 m/s
adat.Vm = adat.Vm.*reshape(repelem(meshgrid(adat.glider_CTD_sos./...
    adat.ADCP_sos,adat.range_cells),1,1,adat.nbeams),adat.ncells,...
    adat.nbeams,length(adat.time));

%% Remove amplitudes> 75dB & SNR<20 (27dB noise floor)
% Following Todd et al 2017
adat.Am(adat.Am>75)=NaN;
adat.Am(10.^(adat.Am-27)<20)=NaN;
adat.Vm(isnan(adat.Am))=NaN;

%% Get rid of correlation values less than 50%
% Following Todd et al 2017
adat.Vm(adat.Cr<50)=NaN;

%% Extract beam information
ds.TrueDepthBeam1=nan(size(squeeze(adat.Vm(:,1,:))));
ds.TrueDepthBeam2=nan(size(squeeze(adat.Vm(:,1,:))));
ds.TrueDepthBeam3=nan(size(squeeze(adat.Vm(:,1,:))));
ds.TrueDepthBeam4=nan(size(squeeze(adat.Vm(:,1,:))));
ds.time = adat.time;
ds.heading = adat.heading_corr;
ds.roll = adat.roll;
ds.pitch = adat.pitch;
ds.depth = adat.glider_depth;
ds.lon   = adat.glider_lon;
ds.lat   = adat.glider_lat;
ds.VelocityBeam1 = squeeze(adat.Vm(:,1,:));
ds.VelocityBeam2 = squeeze(adat.Vm(:,2,:));
ds.VelocityBeam3 = squeeze(adat.Vm(:,3,:));
ds.VelocityBeam4 = squeeze(adat.Vm(:,4,:));
%  Add Bottom Tracking Information
adat.bt_range(abs(adat.bt_vel)>1)=NaN;
adat.bt_vel(abs(adat.bt_vel)>1)=NaN;
ds.BT_time = adat.bt_time;
ds.BT_VelocityBeam1 = squeeze(adat.bt_vel(1,:));
ds.BT_VelocityBeam2 = squeeze(adat.bt_vel(2,:));
ds.BT_VelocityBeam3 = squeeze(adat.bt_vel(3,:));
ds.BT_VelocityBeam4 = squeeze(adat.bt_vel(4,:));
ds.BT_RangeBeam1 = squeeze(adat.bt_range(1,:));
ds.BT_RangeBeam2 = squeeze(adat.bt_range(2,:));
ds.BT_RangeBeam3 = squeeze(adat.bt_range(3,:));
ds.BT_RangeBeam4 = squeeze(adat.bt_range(4,:));

%% Correct range for pitch & roll
% Code from J. Gradone's GITHUB python code
for i = 1:length(adat.time)
    ds.TrueDepthBeam1(:,i) = cell_vert(adat.pitch(i), adat.roll(i),adat.range_cells,1);
    ds.TrueDepthBeam2(:,i) = cell_vert(adat.pitch(i), adat.roll(i),adat.range_cells,2);
    ds.TrueDepthBeam3(:,i) = cell_vert(adat.pitch(i), adat.roll(i),adat.range_cells,3);
    ds.TrueDepthBeam4(:,i) = cell_vert(adat.pitch(i), adat.roll(i),adat.range_cells,4);

    ds.BT_RangeBeam1(i)  = cell_vert(adat.pitch(i), adat.roll(i),ds.BT_RangeBeam1(i),1);
    ds.BT_RangeBeam2(i)  = cell_vert(adat.pitch(i), adat.roll(i),ds.BT_RangeBeam2(i),2);
    ds.BT_RangeBeam3(i)  = cell_vert(adat.pitch(i), adat.roll(i),ds.BT_RangeBeam3(i),3);
    ds.BT_RangeBeam4(i)  = cell_vert(adat.pitch(i), adat.roll(i),ds.BT_RangeBeam4(i),4);
end

%% START LOOP FOR ADCP Processing
% from J. Gradone
dz=2; % depth bin resolution
wDAC = 5; % nr points for weights
wSmoothness = 0.5; % trying 0.5 instead of 1 

% Initialize
out.z_grid = 0:dz:ceil(max(ds.depth,[],'omitnan'));
out.ad2cp_u_dac=NaN(length(adat.gliderSegStartTime),1);
out.ad2cp_v_dac=NaN(length(adat.gliderSegStartTime),1);
out.glider_u_dac=NaN(length(adat.gliderSegStartTime),1);
out.glider_v_dac=NaN(length(adat.gliderSegStartTime),1);
out.glider_lon=NaN(length(adat.gliderSegStartTime),1);
out.glider_lat=NaN(length(adat.gliderSegStartTime),1);
out.seg_mid_time=NaN(length(adat.gliderSegStartTime),1);
out.ad2cp_u_z=NaN(length(out.z_grid),length(adat.gliderSegStartTime));
out.ad2cp_v_z=NaN(length(out.z_grid),length(adat.gliderSegStartTime));

for ii = 1:length(adat.gliderSegStartTime)
    %% Select subset (could be put into loop)
    if (adat.gliderSegEndTime(ii)-adat.gliderSegStartTime(ii))*24*60>10
        ds_subset=get_subset(ds,adat.gliderSegStartTime(ii),adat.gliderSegEndTime(ii));
        if (nanmax(ds_subset.time)-nanmin(ds_subset.time))*24*60>20 %30 s/ping and 20 min dive   

            %% Grid Subset
            % Again code as per Joe Gradone
            ds_subset = grid_ad2cp_beams(ds_subset,adat.range_cells);

            %% Remove bin data that is below glider altimeter depth
            segIdx = adat.time>=min(ds_subset.time) & adat.time<=max(ds_subset.time);
            mean_seg_water_depth = mean(adat.glider_water_depth(segIdx),'omitnan');           
            DepthMask = (ds_subset.TrueDepth+ds_subset.depth)-mean_seg_water_depth>=-1;
            ds_subset.InterpVelocityBeam1(DepthMask)=NaN;
            ds_subset.InterpVelocityBeam2(DepthMask)=NaN;
            ds_subset.InterpVelocityBeam3(DepthMask)=NaN;
            ds_subset.InterpVelocityBeam4(DepthMask)=NaN;

            %% Apply Bottom Tracking Velocities
            % Made this part up. Working on a fix with Clark Richards (DFO-BIO)
            ds_subset.BT_VelocityBeam1(ds_subset.BT_RangeBeam1>15)=NaN;
            ds_subset.BT_VelocityBeam2(ds_subset.BT_RangeBeam2>15)=NaN;
            ds_subset.BT_VelocityBeam3(ds_subset.BT_RangeBeam3>15)=NaN;
            ds_subset.BT_VelocityBeam4(ds_subset.BT_RangeBeam4>15)=NaN;

            res1=(ds_subset.BT_VelocityBeam1-ds_subset.InterpVelocityBeam1);
            res2=(ds_subset.BT_VelocityBeam2-ds_subset.InterpVelocityBeam2);
            res3=(ds_subset.BT_VelocityBeam3-ds_subset.InterpVelocityBeam3);
            res4=(ds_subset.BT_VelocityBeam4-ds_subset.InterpVelocityBeam4);
            ds_subset.InterpVelocityBeam1(~isnan(res1))=...
                ds_subset.InterpVelocityBeam1(~isnan(res1))+res1(~isnan(res1));
            ds_subset.InterpVelocityBeam2(~isnan(res2))=...
                ds_subset.InterpVelocityBeam2(~isnan(res2))+res2(~isnan(res2));
            ds_subset.InterpVelocityBeam3(~isnan(res3))=...
                ds_subset.InterpVelocityBeam3(~isnan(res3))+res3(~isnan(res3));
            ds_subset.InterpVelocityBeam4(~isnan(res4))=...
                ds_subset.InterpVelocityBeam4(~isnan(res4))+res4(~isnan(res4));

            %% Ensemble average (2 ensembles)
            % Skipping this step since it might shift bins in time
            % ds_subset = ensemble_average(ds_subset,2);

            % PLOT Beam 1
%             plotBeam1Profile

            %% Beam to ENU
            % Code from Joe Gradone
            ds_subset=beam2enu(ds_subset,adat.beam2xyz);

            %% Plot Result ENU Transformation
%             plotUVProfile
            
            %% Remove velocities where pitch is less than 10 deg or roll is greater than 5 deg
            % Seems standard
            idxStall = abs(ds_subset.pitch)<10 | abs(ds_subset.roll)>3;
            ds_subset.UVelocity(:,idxStall)=NaN;
            ds_subset.VVelocity(:,idxStall)=NaN;
            ds_subset.WVelocity(:,idxStall)=NaN;

            %% High Velocity Threshold Rel. to STD and Mean
            % Sort of follows Todd et al., 2017 but might be debatable
            %Filter out high velocities relative to glider
            idUV = abs(ds_subset.UVelocity)>abs(mean(ds_subset.UVelocity(:),...
                'omitnan'))+2*abs(std(ds_subset.UVelocity(:),'omitnan'));   
            idVV = abs(ds_subset.VVelocity)>abs(mean(ds_subset.VVelocity(:),...
                'omitnan'))+2*abs(std(ds_subset.VVelocity(:),'omitnan'));
            idWV = abs(ds_subset.WVelocity)>abs(mean(ds_subset.WVelocity(:),...
                'omitnan'))+2*abs(std(ds_subset.WVelocity(:),'omitnan'));
            ds_subset.UVelocity(idUV)=NaN;
            ds_subset.VVelocity(idVV)=NaN;
            ds_subset.WVelocity(idWV)=NaN;
            
            %% PLOT AFTER QC
%             plotUVProfile

            %% Invert ENU to E-W, N-S velocities
            out.glider_u_dac(ii)=mean(adat.glider_u_dac(segIdx),'omitnan');
            out.glider_v_dac(ii)=mean(adat.glider_v_dac(segIdx),'omitnan');
            out.glider_lon(ii)=mean(adat.glider_lon(segIdx),'omitnan');
            out.glider_lat(ii)=mean(adat.glider_lat(segIdx),'omitnan');
            out.seg_mid_time(ii)=mean(ds_subset.time,'omitnan');

            % Inversion from Joe Gradone (Rutgers)
            [O_ls, G_ls, bin_new, obs_per_bin] = ad2cp_inversion(...
                ds_subset.UVelocity, ds_subset.VVelocity, dz, ...
                out.glider_u_dac(ii), out.glider_v_dac(ii), ...
                adat.range_cells, ds_subset.depth, wDAC, wSmoothness);

            out.ad2cp_u_z(:,ii)=interp1(bin_new,real(O_ls),out.z_grid);
            out.ad2cp_v_z(:,ii)=interp1(bin_new,imag(O_ls),out.z_grid);
            out.ad2cp_u_z(abs(out.ad2cp_u_z(:,ii))>0.9,ii)=NaN;
            out.ad2cp_v_z(abs(out.ad2cp_v_z(:,ii))>0.9,ii)=NaN;
            out.ad2cp_u_dac(ii)=mean(out.ad2cp_u_z(:,ii),'omitnan');
            out.ad2cp_v_dac(ii)=mean(out.ad2cp_v_z(:,ii),'omitnan');

%             Plot Resulting Currents
%             figure
%             plot(real(O_ls), bin_new, 'DisplayName', 'u'); hold on
%             plot(imag(O_ls), bin_new, 'DisplayName', 'v');
%             set(gca, 'YDir', 'reverse');
%             title('UV Currents from AD2CP after inversion of ENU')
%             legend('show');
        end
    end
end

save('ADCP_inversion_result.mat','out')

figure(1)
subplot(211)
pcolor(out.seg_mid_time,-out.z_grid,out.ad2cp_u_z);
shading flat
colormap(cmocean('balance'))
cb=colorbar;
caxis([-0.5 0.5])
ylabel(cb,'Glider E-W Currents')
datetick('x','dd.mm','keeplimits')
title('Glider Batray SWOT Calval AD2CP section')


subplot(212)
pcolor(out.seg_mid_time,-out.z_grid,out.ad2cp_v_z);
shading flat
colormap(cmocean('balance'))
cb2=colorbar; 
caxis([-0.5 0.5])
ylabel(cb2,'Glider N-S Currents')
datetick('x','dd.mm','keeplimits')
save_figure(gcf,['./plots/glider_ADCP_section'],[7.5 5],'.png','300')


figure(2)
plot(adat.glider_lon,adat.glider_lat,'.'); hold on;
borders('Canada','facecolor','g');
idnan = ~isnan(adat.glider_u_dac_o);
s1=quiver(adat.glider_lon(idnan),adat.glider_lat(idnan),...
    adat.glider_u_dac(idnan),adat.glider_v_dac(idnan),'filled','Color','r','AutoScale','off','LineWidth',2);
set(gca,'xlim',[-54.9 -50.5],'ylim',[46.7 50.9]);
xlim = get(gca,'XLim'); ylim = get(gca,'Ylim');
title('Glider Batray SWOT Calval Mission')
[x,y,z]=load_bathy(ylim,xlim);
contour(x,y,z','color','k','levellist',[-1000,-300,-200,-100,-50],'showtext','on');
s2=quiver(out.glider_lon,out.glider_lat,out.ad2cp_u_dac,out.ad2cp_v_dac,'filled','Color','m','AutoScale','off');
legend([s1 s2],{'Glider DAC','AD2CP DAC'})
save_figure(gcf,['./plots/glider_ADCP_dac_map'],[7.5 5],'.png','300')
