% Calculates oxygen saturation and concentration

%%



% Calculate oxygen saturation
mooring_data.dox2_sat = mooring_data.dox2_umolkg./(mooring_data.dox2_sol_umolkg.*(mooring_data.atmosphericpress_Pa/constants.atm_in_Pa));

mooring_data.dox2_sat_qc = max(mooring_data.dox2_umolkg_qc,mooring_data.dox2_sol_umolkg_qc);

% Calculate measured oxygen in moles per cubic metre
mooring_data.dox2_molm3 = mooring_data.density_kgm3.*mooring_data.dox2_umolkg.*(10^-6);

mooring_data.dox2_molm3_qc = max(mooring_data.density_kgm3,mooring_data.dox2_umolkg);

% Calculate solubility in moles per cubic metre

mooring_data.dox2_sol_molm3 = mooring_data.density_kgm3.*mooring_data.dox2_sol_umolkg.*(10^-6);

mooring_data.dox2_sol_molm3_qc = max(mooring_data.density_kgm3,mooring_data.dox2_sol_umolkg);