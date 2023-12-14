%% Ecopuck Processing
addpath('./processing_tools/');

% Dark Counts Calc. in-situ
sunfisha_range = sunfish.pressure>350 & sunfish.pressure<900 & ...
    sunfish.time>sunfishetime(2022,01,20) & sunfish.time<sunfishetime(2022,04,01);

% Backscatter 700nm 
eco_bb_dark_counts = 47; % factory value
eco_bb_scale_factor = 1.913E-06; % m-1 sr-1 / counts
% eco_bb_dark_counts = 110; %min_75th_percentile(sunfish.bb700_sig(sunfisha_range))

lambda = 700; % 700 nm backscatter 
theta = 117; % angle of instrument
VBSC  = eco_bb_scale_factor*(sunfish.bb700_sig-eco_bb_dark_counts);  

% 1. Method Zhang et al., 2009 Volume Scattering of Seawater
beta_sw=betasw(lambda, sunfish.temperature, theta, sunfish.salinity_cor);

% Particle backscatter coefficient, BBP 
Xp = 1.1; %recommended by WETLabs 2013
sunfish.bbp700 = 2*pi*(VBSC - beta_sw) * Xp;

% Chlorophyll
eco_chl_dark_counts = 45; % factory value
eco_chl_scale_factor = 0.0073; % micro-g L-1

% Investigate chl_dark_counts from values when readings where dark
eco_chl_dark_counts = min_75th_percentile(sunfish.chlor_sig(sunfisha_range))

sunfish.chlor  = eco_chl_scale_factor*(sunfish.chlor_sig-eco_chl_dark_counts);  

% CDOM
eco_cdom_dark_counts = 50;
eco_cdom_scale_factor = 0.0904;
eco_cdom_dark_counts = min_75th_percentile(sunfish.cdom_sig(sunfisha_range))

sunfish.cdom  = eco_cdom_scale_factor*(sunfish.cdom_sig-eco_cdom_dark_counts);  

% remove weird sunfisha
id = abs(sunfish.roll)>5 | sunfish.chlor<0 | sunfish.bbp700<0 | sunfish.cdom<0;
sunfish.chlor(id)=NaN;
sunfish.bbp700(id)=NaN;
sunfish.cdom(id)=NaN;


