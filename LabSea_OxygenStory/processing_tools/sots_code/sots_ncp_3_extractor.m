clear ; clc; close all

load('./../../l2_corrected_profiles.mat');
load('./../../ERA5_pearldiver.mat');

mooring_data.time = datenum(2020,01,20):1/24:datenum(2020,05,20);
mooring_data.time = mooring_data.time(:);
mooring_data.temp_C = mean_interp(gl.time,nanmean(gl.T(1:10,:),1),mooring_data.time);
mooring_data.psal_PSU = mean_interp(gl.time,nanmean(gl.S(1:10,:),1),mooring_data.time);
mooring_data.dox2_umolkg = mean_interp(gl.time,nanmean(gl.O2_mmolkg(1:10,:),1),mooring_data.time);
mooring_data.mld_m = mean_interp(gl.time,gl.MLD,mooring_data.time);
mooring_data.windspeed_ms = mean_interp(ERA5.time,sqrt(ERA5.u10.^2 + ERA5.v10.^2),mooring_data.time);
mooring_data.atmosphericpress_Pa = mean_interp(ERA5.time,ERA5.mslp,mooring_data.time);
mooring_data.gastension_Pa = flipud(mat2gray(mooring_data.dox2_umolkg))*101325;


mooring_data.temp_C_qc = mooring_data.temp_C;
mooring_data.psal_PSU_qc = mooring_data.psal_PSU;
mooring_data.density_kgm3 = sw_dens(mooring_data.psal_PSU,mooring_data.temp_C,0);
mooring_data.density_kgm3_qc = mooring_data.density_kgm3;
mooring_data.dox2_umolkg_qc = mooring_data.dox2_umolkg;
mooring_data.dox2_sol_umolkg = gsw_O2sol_SP_pt(mooring_data.psal_PSU,mooring_data.temp_C);
mooring_data.dox2_sol_umolkg_qc = mooring_data.dox2_sol_umolkg;
mooring_data.gastension_Pa_qc = mooring_data.gastension_Pa;
mooring_data.mld_m_qc = mooring_data.mld_m;
mooring_data.windspeed_ms_qc = mooring_data.windspeed_ms;
% % 
% % o2_flux_calc(mooring_data,'W09')
% % o2_flux_calc(mooring_data,'S09')
% % o2_flux_calc(mooring_data,'L13')


% function temp_intp=mean_interp(t,v,ti)
%     temp_intp = NaN*ti;
%     t = t(:); v = v(:); ti = ti(:);
%     idx = ~isnan(v);
%     temp_var  = v(idx);
%     temp_x    = t(idx);
%     [~,~,loc]=histcounts(temp_x,ti);
%     temp_x(loc==0)=[]; temp_var(loc==0)=[]; loc(loc==0)=[];           
%     x_mean   = accumarray(loc(:),temp_x(:))  ./accumarray(loc(:),1);
%     var_mean = accumarray(loc(:),temp_var(:))./accumarray(loc(:),1);
%     idnan = find(~isnan(var_mean));
% %     [val,id1]=nanmin(abs(ti-nanmin(x_mean(idnan))));
% %     [val,id2]=nanmin(abs(ti-nanmax(x_mean(idnan))));
%     id1 = loc(1);
%     id2 = loc(end);
%     temp_intp = interp1(x_mean(idnan),var_mean(idnan),ti,'linear','extrap');
% end

% Extracts the data from the relevant netcdf files

% Long term, could all the sots data be compiled into one netcdf file, and
% then selected by a deployment flag? In the meantime, use a cruder
% solution

% function [mooring_data] = sots_ncp_extractor(deployment)
% 
%     if deployment == 'pulse7'
% 
%         % Opens the relevant netcdf file
% 
%         ncid = netcdf.open('IMOS_DWM-SOTS_REFOBKGTPCS_20100817_Pulse_FV02_Pulse-7-2010-Gridded-Data_END-20110708_C-20190522.nc');
% 
%         % In order to access a summary of the information contained in the netcdf
%         % file, you can use the following command
%         % ncdisp('IMOS_DWM-SOTS_REFOBKGTPCS_20100817_Pulse_FV02_Pulse-7-2010-Gridded-Data_END-20110708_C-20190522.nc')
% 
%         % We create a cell object containing the names of each of the variables
%         % stored in the netcdf. Note that variables in a NetCDF file are indexed
%         % from 0.
% 
%         varnames = cell(size(netcdf.inqVarIDs(ncid)));
% 
%         for i=0:(length(netcdf.inqVarIDs(ncid))-1)
% 
%             varnames{i+1} = netcdf.inqVar(ncid,i);
% 
%         end
% 
%         % We create and populate a Matalb structure named "pulse7", containing     
%         % fields named for each of the variables of the NetCDF file, and subfields
%         % containing the instrument name and deployment depth, where relevant.
%         %
%         % For example,          pulse7.DOX2.Aanderaa_Optode3975C_38_5m
% 
%         for i=0:length(varnames)-1
% 
%             try
%                 depths_dummy = netcdf.getAtt(ncid,i,'sensor_depth');
% 
%                 depths_dummy = split(depths_dummy,';');
% 
%                 depths_dummy = strtrim(depths_dummy);
% 
%                 depths_dummy = strcat(depths_dummy,'m');
% 
%                 sensor_dummy = netcdf.getAtt(ncid,i,'sensor_name');
% 
%                 sensor_dummy = split(sensor_dummy,';');
% 
%                 sensor_dummy = strtrim(sensor_dummy);
% 
%                 serial_dummy = netcdf.getAtt(ncid,i,'sensor_serial_number');
% 
%                 serial_dummy = split(serial_dummy,';');
% 
%                 serial_dummy = strtrim(serial_dummy);
% 
%                 names_dummy = matlab.lang.makeValidName(strcat(sensor_dummy,'_',serial_dummy,'_',depths_dummy));
% 
%                 dummy_variable = netcdf.getVar(ncid,i);
% 
%                 for j=1:length(names_dummy)
% 
%                     pulse7.(varnames{i+1}).(names_dummy{j}) = dummy_variable(:,j);
% 
%                 end
% 
% 
%             catch
% 
%                 pulse7.(varnames{i+1}) = netcdf.getVar(ncid,i);
% 
%             end
% 
%         end
% 
%         % We transpose the time vector, as it has a different format.
% 
%         pulse7.TIME = pulse7.TIME;
% 
%         % Here we set a time offset, as MATLAB's serial date number format is based
%         % on time since year 0, while the netcdf data's serial date number is from
%         % 01/01/1950 00:00:00 UTC
% 
%         time_offset = datenum(1950,1,1,0,0,0);
% 
%         % We apply this time offset to the TIME variable of the structure, so that
%         % data is now correctly timestamped in Matlab.
% 
%         pulse7.TIME = pulse7.TIME + time_offset.*ones(size(pulse7.TIME));
% 
%         % We close the netcdf file
% 
%         netcdf.close(ncid);
% 
%         % We now assign the relevant data from the extraction to the outputs of
%         % the function
% 
%         mooring_data.time = pulse7.TIME;
% 
%         mooring_data.temp_C = pulse7.TEMP.Sea_BirdElectronics_SBE16plusV2_6331_31_1m;
% 
%         mooring_data.temp_C_qc = pulse7.TEMP_quality_control(:,1);
% 
%         mooring_data.psal_PSU = pulse7.PSAL.Sea_BirdElectronics_SBE16plusV2_6331_31_1m;
% 
%         mooring_data.psal_PSU_qc = pulse7.PSAL_quality_control(:,1);
% 
%         mooring_data.density_kgm3 = pulse7.DENSITY.Sea_BirdElectronics_SBE16plusV2_6331_31_1m;
% 
%         mooring_data.density_kgm3_qc = pulse7.DENSITY_quality_control(:,1);
% 
%         % DOX2.Aanderaa_Optode3975_1161_31_1m or Sea_BirdElectronics_SBE43_431635_31_1m
%         mooring_data.dox2_umolkg = pulse7.DOX2.Aanderaa_Optode3975_1161_31_1m; 
% 
%         mooring_data.dox2_umolkg_qc = pulse7.DOX2_quality_control(:,1);
% 
%         mooring_data.dox2_sol_umolkg = pulse7.OXSOL.Sea_BirdElectronics_SBE16plusV2_6331_31_1m;
% 
%         mooring_data.dox2_sol_umolkg_qc = pulse7.OXSOL_quality_control(:,1);
% 
%         mooring_data.gastension_Pa = pulse7.TOTAL_GAS_PRESSURE.ProOceanus_GTD_29_102_15_31_1m*100; % Record is in hPa
% 
%         mooring_data.gastension_Pa_qc = pulse7.TOTAL_GAS_PRESSURE_quality_control;
% 
%         mooring_data.mld_m = pulse7.MLD;
% 
%         mooring_data.mld_m_qc = pulse7.MLD_quality_control;
% 
%         % Use Eric's simulated wind speeds (Eric Schulz BOM 10m winds)
%         load('Pulse_7_Eric_wind_sim.mat','wind_sim_total');
% 
%         mooring_data.windspeed_ms = wind_sim_total;
% 
%     
% 
%     elseif deployment == 'pulse9'
% 
%         % Here the NetCDF file containing the Pulse 9 data is imported into
%         % Matlab. This file must be in Matlab's current folder.
% 
%         ncid = netcdf.open('IMOS_ABOS-SOTS_RWBKGOTPCS_20120619_Pulse_FV02_Pulse-9-2012-Gridded-Data_END-20130510_C-20210131.nc');
% 
% 
%         % In order to access a summary of the information contained in the netcdf
%         % file, you can use the following command
%         % ncdisp('IMOS_ABOS-SOTS_RWBKGOTPCS_20120619_Pulse_FV02_Pulse-9-2012-Gridded-Data_END-20130510_C-20190910.nc')
% 
%         % We create a cell object containing the names of each of the variables
%         % stored in the netcdf. Note that variables in a NetCDF file are indexed
%         % from 0.
% 
%         varnames = cell(size(netcdf.inqVarIDs(ncid)));
% 
%         for i=0:(length(netcdf.inqVarIDs(ncid))-1)
% 
%             varnames{i+1} = netcdf.inqVar(ncid,i);
% 
%         end
% 
%         % We create and populate a Matalb structure named "pulse9", containing     
%         % fields named for each of the variables of the NetCDF file, and subfields
%         % containing the instrument name and deployment depth, where relevant.
%         %
%         % For example,          pulse9.DOX2.Aanderaa_Optode3975C_38_5m
% 
%         for i=0:length(varnames)-1
% 
%             try
%                 depths_dummy = netcdf.getAtt(ncid,i,'sensor_depth');
% 
%                 depths_dummy = split(depths_dummy,';');
% 
%                 depths_dummy = strtrim(depths_dummy);
% 
%                 depths_dummy = strcat(depths_dummy,'m');
% 
%                 sensor_dummy = netcdf.getAtt(ncid,i,'sensor_name');
% 
%                 sensor_dummy = split(sensor_dummy,';');
% 
%                 sensor_dummy = strtrim(sensor_dummy);
% 
%                 names_dummy = matlab.lang.makeValidName(strcat(sensor_dummy,'_',depths_dummy));
% 
%                 dummy_variable = netcdf.getVar(ncid,i);
% 
%                 for j=1:length(names_dummy)
% 
%                     pulse9.(varnames{i+1}).(names_dummy{j}) = dummy_variable(:,j);
% 
%                 end
% 
% 
%             catch
% 
%                 dummy_variable = netcdf.getVar(ncid,i);
% 
%                 if length(dummy_variable) > 100
% 
%                     pulse9.(varnames{i+1}) = dummy_variable; 
% 
%                 else
% 
%                     pulse9.(varnames{i+1}) = netcdf.getVar(ncid,i);
% 
%                 end
% 
%             end
% 
%         end
% 
%         % We transpose the time vector, as it has a different format.
% 
%         pulse9.TIME = pulse9.TIME;
% 
%         % Here we set a time offset, as MATLAB's serial date number format is based
%         % on time since year 0, while the netcdf data's serial date number is from
%         % 01/01/1950 00:00:00 UTC
% 
%         time_offset = datenum(1950,1,1,0,0,0);
% 
%         % We apply this time offset to the TIME variable of the structure, so that
%         % data is now correctly timestamped in Matlab.
% 
%         pulse9.TIME = pulse9.TIME + time_offset.*ones(size(pulse9.TIME));
% 
%         % We adjust the VEMCO temperature sensors from 45-85m in order for them to
%         % read in agreement with the two Seabird sensors they sit between, by
%         % applying linear offsets.
% 
%         vemco_offsets = [0.025 0.025 0.015 0.045 0.015 -0.02 0.03 0];
% 
%         fields_dummy=fieldnames(pulse9.TEMP);
% 
%         for i=1:length(vemco_offsets)
% 
%             pulse9.TEMP.(fields_dummy{i+1}) = pulse9.TEMP.(fields_dummy{i+1}) + vemco_offsets(i);
% 
%         end
% 
%         % We close the netcdf file
% 
%         netcdf.close(ncid);
% 
%         % We import the matching sofs3 surface data into Matlab, to access
%         % windspeed data.
% 
%         ncid2 = netcdf.open('IMOS_DWM-ASFS_CFMST_20120714T080000Z_SOFS_FV02_SOFS-3-2012_END-20130102T000000Z_C-20180716T064518Z.nc');
% 
%         % We perform the same routine as with the Pulse 9 data, to extract the
%         % data and place it in a structure named "sofs3".
% 
%         varnames2 = cell(size(netcdf.inqVarIDs(ncid2)));
% 
%         for i=0:(length(netcdf.inqVarIDs(ncid2))-1)
% 
%             varnames2{i+1} = netcdf.inqVar(ncid2,i);
% 
%         end
% 
%         for i=0:length(varnames2)-1
% 
%             try
%                 depths_dummy = netcdf.getAtt(ncid2,i,'sensor_depth');
% 
%                 depths_dummy = split(depths_dummy,';');
% 
%                 depths_dummy = strtrim(depths_dummy);
% 
%                 depths_dummy = strcat(depths_dummy,'m');
% 
%                 sensor_dummy = netcdf.getAtt(ncid2,i,'sensor_name');
% 
%                 sensor_dummy = split(sensor_dummy,';');
% 
%                 sensor_dummy = strtrim(sensor_dummy);
% 
%                 names_dummy = matlab.lang.makeValidName(strcat(sensor_dummy,'_',depths_dummy));
% 
%                 dummy_variable = netcdf.getVar(ncid2,i);
% 
%                 for j=1:length(names_dummy)
% 
%                     sofs3.(varnames2{i+1}).(names_dummy{j}) = dummy_variable(:,j);
% 
%                 end
% 
% 
%             catch
% 
%                 sofs3.(varnames2{i+1}) = netcdf.getVar(ncid2,i);
% 
%             end
% 
%         end
% 
%         sofs3.TIME = sofs3.TIME + time_offset.*ones(size(sofs3.TIME));
% 
%         % As sofs3 surface data has been recorded in minute intervals, from a
%         % different starting date to Pulse 10, here we create a time index to
%         % subsample the sofs3 data to align it with Pulse 9 sampling.
% 
%         % This results in the following being true:
%         % sofs3.TIME(sofs3timeindex) = pulse9.TIME
% 
%         sofs3timeindex = zeros(size(pulse9.TIME));
% 
%         for i = 1:length(pulse9.TIME)
% 
%             try
% 
%                 sofs3timeindex(i) = find(sofs3.TIME == pulse9.TIME(i),1);
% 
%             catch
% 
%                 sofs3timeindex(i) = 1;
% 
%             end
% 
%         end
% 
%         % We close the netcdf
%         netcdf.close(ncid2);
% 
%         % We now assign the relevant data from the extraction to the outputs of
%         % the function
% 
%         mooring_data.time = pulse9.TIME;
% 
%         mooring_data.temp_C = pulse9.TEMP.Sea_BirdElectronics_SBE16plusV2_28_5m;
% 
%         mooring_data.temp_C_qc = pulse9.TEMP_quality_control(:,1);
% 
%         mooring_data.psal_PSU = pulse9.PSAL.Sea_BirdElectronics_SBE16plusV2_28_5m;
% 
%         mooring_data.psal_PSU_qc = pulse9.PSAL_quality_control(:,1);
% 
%         mooring_data.density_kgm3 = pulse9.DENSITY.Sea_BirdElectronics_SBE16plusV2_28_5m;
% 
%         mooring_data.density_kgm3_qc = pulse9.DENSITY_quality_control(:,1);
% 
%         % Aanderaa_Optode3975C_38_5m or Sea_BirdElectronics_SBE43_38_5m
%         mooring_data.dox2_umolkg = pulse9.DOX2.Aanderaa_Optode3975C_28_5m; 
% 
%         mooring_data.dox2_umolkg_qc = pulse9.DOX2_quality_control(:,1);
% 
%         mooring_data.dox2_sol_umolkg = pulse9.OXSOL.Sea_BirdElectronics_SBE16plusV2_28_5m;
% 
%         mooring_data.dox2_sol_umolkg_qc = pulse9.OXSOL_quality_control(:,1);
% 
%         mooring_data.gastension_Pa = pulse9.TOTAL_GAS_PRESSURE.ProOceanus_GTD_28_5m/(1E3); % Record is in millibars in netcdf
% 
%         mooring_data.gastension_Pa_qc = pulse9.TOTAL_GAS_PRESSURE_quality_control;
% 
%         mooring_data.mld_m = pulse9.MLD;
% 
%         mooring_data.mld_m_qc = pulse9.MLD_quality_control;
% 
%         mooring_data.windspeed_ms = sofs3.WSPD10M(sofs3timeindex);
% 
%         % No atmospheric pressure available
% 
%     
% 
%     elseif deployment == 'sofs75'
% 
% 
%         %     addpath('C:\xx\xx\data_folder\composite1.dat')
%         %     
%         %     addpath('C:\xx\xx\data_folder\composite1.dat')
% 
%         % Here the NetCDF file containing the SOFS-7.5 data is imported into
%         % Matlab. This file must be in Matlab's current folder.
% 
%         ncid = netcdf.open('IMOS_ABOS-SOTS_FRVMWEBUOSTCGP_20180628_SOFS_FV02_SOFS-7.5-2018-Gridded-Data_END-20190423_C-20190916.nc');
% 
%         % In order to access a summary of the information contained in the netcdf
%         % file, you can use the following command
%         % ncdisp('IMOS_ABOS-SOTS_FRVMWEBUOSTCGP_20180628_SOFS_FV02_SOFS-7.5-2018-Gridded-Data_END-20190423_C-20190916.nc')
% 
%         % We create a cell object containing the names of each of the variables
%         % stored in the netcdf. Note that variables in a NetCDF file are indexed
%         % from 0.
% 
%         varnames = cell(size(netcdf.inqVarIDs(ncid)));
% 
%         for i=0:(length(netcdf.inqVarIDs(ncid))-1)
% 
%             varnames{i+1} = netcdf.inqVar(ncid,i);
% 
%         end
% 
%         % We create and populate a Matalb structure named "sofs75", containing     
%         % fields named for each of the variables of the NetCDF file, and subfields
%         % containing the instrument name and deployment depth, where relevant.
%         %
%         % For example,          sofs75.DOX2.Aanderaa_Optode3975C_38_5m
% 
%         for i=0:length(varnames)-1
% 
%             try
%                 depths_dummy = netcdf.getAtt(ncid,i,'sensor_depth');
% 
%                 depths_dummy = split(depths_dummy,';');
% 
%                 depths_dummy = strtrim(depths_dummy);
% 
%                 depths_dummy = strcat(depths_dummy,'m');
% 
%                 sensor_dummy = netcdf.getAtt(ncid,i,'sensor_name');
% 
%                 sensor_dummy = split(sensor_dummy,';');
% 
%                 sensor_dummy = strtrim(sensor_dummy);
% 
%                 serial_dummy = netcdf.getAtt(ncid,i,'sensor_serial_number');
% 
%                 serial_dummy = split(serial_dummy,';');
% 
%                 serial_dummy = strtrim(serial_dummy);
% 
%                 names_dummy = matlab.lang.makeValidName(strcat(sensor_dummy,'_',serial_dummy,'_',depths_dummy));
% 
%                 dummy_variable = netcdf.getVar(ncid,i);
% 
%                 for j=1:length(names_dummy)
% 
%                     sofs75.(varnames{i+1}).(names_dummy{j}) = dummy_variable(:,j);
% 
%                 end
% 
% 
%             catch
% 
%                 sofs75.(varnames{i+1}) = netcdf.getVar(ncid,i);
% 
%             end
% 
%         end
% 
%         % We extract the global attributes from the netcdf file and insert them
%         % into the structure
% 
% 
% 
%         % We transpose the time vector, as it has a different format.
% 
%         sofs75.TIME = sofs75.TIME;
% 
%         % Here we set a time offset, as MATLAB's serial date number format is based
%         % on time since year 0, while the netcdf data's serial date number is from
%         % 01/01/1950 00:00:00 UTC
% 
%         time_offset = datenum(1950,1,1,0,0,0);
% 
%         % We apply this time offset to the TIME variable of the structure, so that
%         % data is now correctly timestamped in Matlab.
% 
%         sofs75.TIME = sofs75.TIME + time_offset.*ones(size(sofs75.TIME));
% 
%         % We close the netcdf file
% 
%         netcdf.close(ncid);
% 
% 
%         % We now assign the relevant data from the extraction to the outputs of
%         % the function
% 
%         mooring_data.time = sofs75.TIME;
% 
%         mooring_data.temp_C = sofs75.TEMP.Sea_BirdElectronics_SBE37SMP_ODO_15696_30_0m;
% 
%         mooring_data.temp_C_qc = sofs75.TEMP_quality_control(:,6);
% 
%         mooring_data.psal_PSU = sofs75.PSAL.Sea_BirdElectronics_SBE37SMP_ODO_15696_30_0m;
% 
%         mooring_data.psal_PSU_qc = sofs75.PSAL_quality_control(:,4);
% 
%         % Need to manually calculate density!
% 
%         SP = mooring_data.psal_PSU;
% 
%         t = mooring_data.temp_C;
% 
%         p = sofs75.PRES.Sea_BirdElectronics_SBE37SMP_ODO_15696_30_0m;
% 
%         long = sofs75.LONGITUDE;
% 
%         lat = sofs75.LATITUDE;
% 
%         SA = gsw_SA_from_SP(SP,p,long,lat);
% 
%         CT = gsw_CT_from_t(SA,t,p);
% 
%         mooring_data.density_kgm3 = gsw_rho(SA,CT,p);
% 
%         mooring_data.density_kgm3_qc = max([mooring_data.psal_PSU_qc,mooring_data.temp_C_qc,sofs75.PRES_quality_control(:,1)],[],2);
% 
%         mooring_data.dox2_umolkg = sofs75.DOX2.Sea_BirdElectronics_SBE37SMP_ODO_15696_30_0m; 
% 
%         mooring_data.dox2_umolkg_qc = sofs75.DOX2_quality_control(:,2);
% 
%         mooring_data.dox2_sol_umolkg = sofs75.OXSOL.Sea_BirdElectronics_SBE37SMP_ODO_15696_30_0m;
% 
%         mooring_data.dox2_sol_umolkg_qc = sofs75.OXSOL_quality_control(:,2);
%         
%         mooring_data.sub_mld_dox2_umolkg = sofs75.DOX2.Sea_BirdElectronics_SBE37SMP_ODO_15972_480_0m;
%         
%         mooring_data.sub_mld_dox2_umolkg_qc = sofs75.DOX2_quality_control(:,end);
% 
%         mooring_data.gastension_Pa = sofs75.TOTAL_GAS_PRESSURE.ProOceanus_TGTD_37_468_33_30_0m*100; % Record is in hPa in netcdf
% 
%         mooring_data.gastension_Pa_qc = sofs75.TOTAL_GAS_PRESSURE_quality_control(:,2);
% 
%         % Build mld from FV02 file
%         build_mld = true;
% 
%         % Calculate sensor offsets from FV02 file
%         manual_offset = true;
% 
%         if build_mld
% 
%             if manual_offset
% 
%                 % If temperature offsets haven't been applied in the FV02 file
%                 % already
%                 sofs75.TEMParray = struct2array(sofs75.TEMP);
% 
%                 % We select the profiles where all instruments record QC values of 1
% 
%                 sofs75.valid_record_start = find(all(sofs75.TEMP_quality_control==1,2),1);
% 
%                 sofs75.valid_record_end = find(all(sofs75.TEMP_quality_control==1,2),1,'last');
% 
%                 sofs75.TEMParray_qc1_profiles = sofs75.TEMParray(sofs75.valid_record_start:sofs75.valid_record_end,:);
% 
%                 sofs75.TIME_qc1_profiles = sofs75.TIME(sofs75.valid_record_start:sofs75.valid_record_end);
% 
%                 % We create an array of Nan the same size as 'sofs75.TEMParray_qc1_profiles', to
%                 % store the excursions of the Starmon sensor values from the Seabird sensor
%                 % values.
% 
%                 sofs75.TEMPexcursions = NaN(size(sofs75.TEMParray_qc1_profiles));
% 
%                 % We set a temperature threshold of 0.004°C, based on the accuracy of
%                 % +-0.002°C of the SBE 37-SMP-ODO, as in a body of water of one temperature,
%                 % two SBE 37s could read 0.004°C out from each other.
% 
%                 temp_thresh = 0.004;
% 
%                 % We give the indices of the trusted temperature sources (generally SBEs)
% 
%                 sofs75.trusted_temps = [1 6 14 18 25];
% 
%                 % We loop through all the selected temperature profiles
% 
%                 for i = 1:length(sofs75.TEMPexcursions)
% 
%                         % We create multiple if conditions. Each condition checks for the current 
%                         % time stamp if the pairs of Seabird sensors closest to each other 
%                         % agree to within the temp_thresh. If they do, the excursions of
%                         % the Starmon sensors in between the two Seabirds are recorded in
%                         % the array 'sofs75.TEMPexcursions' for the relevant time stamp.
% 
%                     for j=1:(length(sofs75.trusted_temps)-1)
% 
%                         if abs(sofs75.TEMParray_qc1_profiles(i,sofs75.trusted_temps(j))-sofs75.TEMParray_qc1_profiles(i,sofs75.trusted_temps(j+1)))<=temp_thresh
% 
%                             sofs75.TEMPexcursions(i,(sofs75.trusted_temps(j)+1):(sofs75.trusted_temps(j+1)-1)) = mean([sofs75.TEMParray_qc1_profiles(i,sofs75.trusted_temps(j)) sofs75.TEMParray_qc1_profiles(i,sofs75.trusted_temps(j+1))])  - sofs75.TEMParray_qc1_profiles(i,(sofs75.trusted_temps(j)+1):(sofs75.trusted_temps(j+1)-1));
% 
%                         end
% 
%                     end
% 
%                 end
%                 % We calculate the median value of the offsets for each temperature sensor,
%                 % ignoring any NaN values generated by the Seabird sensors.
% 
%                 sofs75.TEMPoffsets = median(sofs75.TEMPexcursions,'omitnan');
% 
%                 % We set the NaN results from the Seabirds to 0.
% 
%                 sofs75.TEMPoffsets(isnan(sofs75.TEMPoffsets))=0;
% 
%                 % We adjust the selected temperature readings in sofs75.TEMParray_qc1_profiles 
%                 % and sofs75.TEMParray by the calculated offsets.
% 
%                 sofs75.TEMParray_qc1_profiles = sofs75.TEMParray_qc1_profiles + sofs75.TEMPoffsets;
% 
%                 sofs75.TEMParray = sofs75.TEMParray + sofs75.TEMPoffsets;
% 
%             else
% 
%                 sofs75.TEMParray = struct2array(sofs75.TEMP);
% 
%             end
% 
%             % *********************************************************************** %
%             % ----------------------------------------------------------------------- %
%             %                Calculation of Mixed layer depth                         %
%             % ----------------------------------------------------------------------- %
%             % *********************************************************************** %
% 
%             % We calculate the mixed layer depth using a threshold method. 
% 
%             % First a temperature threshold is set.
% 
%             MLDP_temp_thresh = 0.2;
% 
%             % We create a new field in the sofs75 structure to store the MLDP data.
% 
%             sofs75.MLDP.thresh02 = max(sofs75.DEPTH_TEMP).*ones(length(sofs75.TEMParray),1);
% 
%             % For each timestamp, the MLDP is set at the depth of the first sensor whose
%             % absolute difference from the shallowest sensor exceeds the threshold. If
%             % no sensor exceeds this threshold, then the MLDP is set as the deepest
%             % sensor depth.
% 
%             for i = 1:length(sofs75.MLDP.thresh02)
% 
%                 if any(abs(sofs75.TEMParray(i,1)-sofs75.TEMParray(i,:))>=MLDP_temp_thresh)
% 
%                     sofs75.MLDP.thresh02(i) = sofs75.DEPTH_TEMP(find(abs(sofs75.TEMParray(i,1)-sofs75.TEMParray(i,:))>=MLDP_temp_thresh,1));
% 
%                 end
% 
%             end
% 
%             % Here we use a linear interpolation method, as suggested in Huang et al (2018)
%             % as the most accurate method for temperature data, spare in depth resolution.
% 
%             sofs75.MLDP.interp02 = nan.*ones(length(sofs75.TEMParray),1);
% 
%             MLDP_temp_interp02 = 0.2;
% 
%             for i = 1:length(sofs75.MLDP.interp02)
% 
%                 if sum(isnan(sofs75.TEMParray(i,:)))==0 && ~any(abs(sofs75.TEMParray(i,1)-sofs75.TEMParray(i,:))>=MLDP_temp_interp02)
% 
%                     sofs75.MLDP.interp02(i) = max(sofs75.DEPTH_TEMP);   
% 
%                 elseif sum(isnan(sofs75.TEMParray(i,:)))==0
% 
%                     idx_below = find(abs(sofs75.TEMParray(i,1)-sofs75.TEMParray(i,:))>=MLDP_temp_interp02,1);
% 
%                     idx_above = idx_below-1;
% 
%                     sofs75.MLDP.interp02(i) = interp1(sofs75.TEMParray(i,idx_above:idx_below),sofs75.DEPTH_TEMP(idx_above:idx_below),sofs75.TEMParray(i,1)-MLDP_temp_interp02);        
% 
%                 end
% 
%             end
% 
% 
%         end
% 
%         mooring_data.mld_m = sofs75.MLDP.interp02;
% 
%         mooring_data.mld_m_qc = max(sofs75.TEMP_quality_control,[],2);
% 
%         % Windspeed measured at 2.61m, scale using the power law from
%         % Justus and Mikhail (1976)
%         
%         for w_idx=1:length(sofs75.WSPD(:,1))
%             
%             wind_alpha(w_idx) = (0.37-0.088*log(sofs75.WSPD(w_idx,1)))/(1-0.088*log(2.61/10));
%             
%             windspeed_10m(w_idx) = sofs75.WSPD(w_idx,1) * (10/2.61)^wind_alpha(w_idx);
%             
%         end
%         
% %         % Currently using linearly interpolated 10m windspeeds from the 2
% %         % hour realtime SOTS files, downloaded from the AODN as of
% %         % 14/05/20, will progress to using the 10m windspeeds from FV02
% %         % file once available
% %         load('sofs75_wspd10_interp.mat','WSPD10M_interp')
% %         
% %         mooring_data.windspeed_ms = WSPD10M_interp;
% 
%         mooring_data.windspeed_ms = sofs75.WSPD(:,1);
% 
%         mooring_data.windspeed_ms_qc = sofs75.WSPD_quality_control(:,1);
% 
%         mooring_data.atmosphericpress_Pa = sofs75.CAPH(:,1)*100; % Record is in hPa
% 
%     end
%     
% end


