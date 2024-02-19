% Loads the constants used throughout the SOTS NCP calculation
% Constants are listed under the script in which they are first used
%% sots_ncp_6_O2_converter
constants.atm_in_Pa = 101325;

%% sots_ncp_8_gas_records

% N2 constants from Table 4 of Hamme and Emerson 2004
constants.A0_n2 = 6.42931;
constants.A1_n2 = 2.92704;
constants.A2_n2 = 4.32531;
constants.A3_n2 = 4.69149;
constants.B0_n2 = -7.44129e-3;
constants.B1_n2 = -8.02566e-3;
constants.B2_n2 = -1.46775e-2;

% Wanger polynomial coefficients from Wagner and Pruss (2002)
constants.wagner_a1 = -7.85951783;
constants.wagner_a2 = 1.84408259;
constants.wagner_a3 = -11.7866497;
constants.wagner_a4 = 22.6807411;
constants.wagner_a5 = -15.9618719;
constants.wagner_a6 = 1.80122502;

% Millero (1974) coefficients in Dickson, Sabine et al. (2007)
constants.millero_c0 = 0.90799;
constants.millero_c1 = -0.08992;
constants.millero_c2 = 0.18458;
constants.millero_c3 = -0.07395;
constants.millero_c4 = -0.00221;

% Atmospheric gas fractions from Gluekauf (1951)
constants.mole_fraction_N2 = 0.78084;
constants.mole_fraction_Ar = 0.00934;
constants.mole_fraction_CO2 = 0.000387;
constants.mole_fraction_O2 = 0.20946;

%% sots_ncp_9_gas_transfer_velocity_X

% Schmidt constants from Wanninkhof (14)
%Oxygen
constants.Sc_A_O2_W14 = 1920.4;
constants.Sc_B_O2_W14 = -135.6;
constants.Sc_C_O2_W14 = 5.2122;
constants.Sc_D_O2_W14 = -0.10939;
constants.Sc_E_O2_W14 = 0.00093777;

%Nitrogen
constants.Sc_A_N2_W14 = 2304.8;
constants.Sc_B_N2_W14 = -162.75;
constants.Sc_C_N2_W14 = 6.2557;
constants.Sc_D_N2_W14 = -0.13129;
constants.Sc_E_N2_W14 = 0.0011255;

% Schmidt constants from Wanninkhof (92)
%Oxygen
constants.Sc_A_O2_W92 = 1953.4;
constants.Sc_B_O2_W92 = 128.00;
constants.Sc_C_O2_W92 = 3.9918;
constants.Sc_D_O2_W92 = 0.050091;

%Nitrogen
constants.Sc_A_N2_W92 = 2206.1;
constants.Sc_B_N2_W92 = 144.86;
constants.Sc_C_N2_W92 = 4.5413;
constants.Sc_D_N2_W92 = 0.056988;

% Schmidt number for conversion
constants.schmidt_CO2_20C_35PSU = 660;

%% sots_ncp_10_N2_excahnge

% Molecular diffusion, using Ferrell and Himmelblau (1967).
% Activation energies
constants.Ea_O2 = -18.70; 
constants.Ea_N2 = -18.50;

% Maximum diffusion coefficients
constants.Ac_O2 = 4286;
constants.Ac_N2 = 3412;

% Gas constant
constants.gas_constant = 8.3144621;

% Salinity adjustment
constants.Jahne_sal_co = 0.049;

%% sots_ncp_13_net_com_prod_calc
% This is the NCP ratio of carbon to O2 from Anderson and Sarmiento (1994)
constants.ncp_oxygen2carbon = 1.45; %Redfield AS 1994

constants.atomic_mass_C = 12.0107;