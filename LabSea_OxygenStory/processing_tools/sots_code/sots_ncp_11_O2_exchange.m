 
%% O2 exchange and biology calculation
% Changes in measured O2 between timestamps are examined, using modelling 
% to estimate the amounts of change that can be attributed to physical processes
% and therefore attributing the remainder of the change to biological processes.

% Sub mld oxygen is selected according to the user's choice in the driver code

if ~isfield(mooring_data,'sub_mld_dox2_umolkg') || sub_mld_dox2_manual_override
    
    mooring_data.sub_mld_dox2_molm3 = sub_mld_dox2_choice * ones(size(mooring_data.time));
    
    disp(['Sub MLD O2: user specified constant of ',num2str(sub_mld_dox2_choice),'mol/m^3 used.'])
    
else 
    
    mooring_data.sub_mld_dox2_molm3 = mooring_data.sub_mld_dox2_umolkg.*mooring_data.density_kgm3/(1E6);
   
    disp('Sub MLD O2: Timeseries available from mooring used')
    
end
 
% Zero arrays are preallocated for use in the for loop
 
[mooring_data.GE_dox2_molm2, mooring_data.dox2_gas_exchange_molm2, mooring_data.dox2_bubbles_molm2] = deal(zeros(size(mooring_data.time)));
 

% A zero array representing biological contributions in mols per square metre is created

mooring_data.dox2_biology_molm2 = zeros(size(mooring_data.dox2_molm3));
 
% We loop through the timeseries, estimating the contribution that each
% physical processes will have made to the change in oxygen over the previous
% hour, and finally assign any remaining change not explained by these physical
% contributions to biology. This is a rearrangement of eqn 10, as described 
% in Emerson et al. 2008

for i = 2:length(mooring_data.time)
    
    
    % The rate of gas exchange in mol per square metre per second is calculated
    
    mooring_data.GE_dox2_molm2(i-1) = (mooring_data.Gc_O2_ms(i-1).*(mooring_data.dox2_molm3(i-1) - (mooring_data.atmosphericpress_Pa(i-1)/constants.atm_in_Pa)*mooring_data.dox2_sol_molm3(i-1)));
    
    % The effect of the gas exchange is calculated at the current time stamp 
    % (the multiplication by 3600 here account for the fact that the rates 
    % are per second, and the time stamps are an hour long).
    
    mooring_data.dox2_gas_exchange_molm2(i) = 3600*(-mooring_data.GE_dox2_molm2(i-1));
    
    % The effect of bubble injection at the current time stamp is calculated, 
    % using the calculated Vinj and user defined bubble_beta.
    
    mooring_data.dox2_bubbles_molm2(i) = 3600*mooring_data.Vinj(i-1)*constants.mole_fraction_O2*(1+1/bubble_beta);
    
    % Biological oxygen is calculated by subtracting the calculated physical changes
    % from the measured change, and assigning the remaining change to biology, 
    % as described in Emerson et al. 2008 eqn 10.
    
    mooring_data.dox2_biology_molm2(i) = mooring_data.mld_smooth(i)*(mooring_data.dox2_molm3(i)-mooring_data.dox2_molm3(i-1)) - mooring_data.dox2_gas_exchange_molm2(i) - mooring_data.dox2_bubbles_molm2(i);
    
end
 
%% Physical oxygen model
 % Here we calculate a theoretical 'physical' oxygen record - beginning with
% the first timestamp's measured O2 value, and calculating how we believe
% the O2 would change over time, without any biology being present. 

% We preallocate zero arrays used in the for loop

[mooring_data.GE_dox2_phys_molm3, mooring_data.dox2_phys_gas_exchange_molm3, mooring_data.dox2_phys_bubbles_molm3, mooring_data.dox2_phys_molm3] = deal(zeros(size(mooring_data.time)));

% We set the first value of dox2_phys_molm3 to the first measured oxygen value.

mooring_data.dox2_phys_molm3(1) = mooring_data.dox2_molm3(1);

for i = 2:length(mooring_data.time)

    % We calculate the rate of gas exchange in the previous time stamp, as
    % seen in the final equations of section 8, but now using the
    % theoretical physical oxygen record.
    
    mooring_data.GE_dox2_phys_molm2(i-1) = (mooring_data.Gc_O2_ms(i-1)*((mooring_data.Sc_O2(i-1)/600).^(-0.5)).*(mooring_data.dox2_phys_molm3(i-1) - (mooring_data.atmosphericpress_Pa(i-1)/constants.atm_in_Pa)*mooring_data.dox2_sol_molm3(i-1)));

    

    % We calculate the effect of the above calculated gas exchange in the
    % current time stamp (the multiplication by 3600 here account for the
    % fact that the rates are per second, and the time stamps are an hour
    % long).

    mooring_data.dox2_phys_gas_exchange_molm3(i) = (3600/mooring_data.mld_smooth(i))*(-mooring_data.GE_dox2_phys_molm2(i-1));

    % We calculate the effect of bubble injection in the current time
    % stamp, using the calculated Vinj and bubble_beta.

    mooring_data.dox2_phys_bubbles_molm3(i) = 3600/mooring_data.mld_smooth(i)*mooring_data.Vinj(i-1)*constants.mole_fraction_O2*(1+1/bubble_beta);

    % For each time stamp, we add the calculated contributions to calculate
    % the theoretical physical oxygen record, without entrainment.

    mooring_data.dox2_phys_molm3(i) = mooring_data.dox2_phys_molm3(i-1) + mooring_data.dox2_phys_bubbles_molm3(i) + mooring_data.dox2_phys_gas_exchange_molm3(i);

end

% We convert the timeseries from mol per cubic metre to micromoles per kilogram

mooring_data.dox2_phys_umkg = (mooring_data.dox2_phys_molm3*1E6)./(mooring_data.density_kgm3);








