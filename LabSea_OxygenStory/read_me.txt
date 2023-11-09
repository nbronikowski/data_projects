%% Folder Contains Steps to Process / Extract Figures for Oxygen Paper

% Part 1: Extract Glider Data from ERDDAP

% Part 2: Data Preparation and Cleaning
% - basic qc
% - clean salinity and correct for thermal lag
% - data gridding
% - correct oxygen for response time lag
% - compare T,S, O2 to nearby mooring/ship data to apply a calibration
% - calculate optode drift from near-surface measurements

% Part 3: Compute additional information
% - mixed layer depth
% - chlorophyll, backscatter, cdom for Sunfish only
% - potentially apply quenching correction for Sunfish data

% Part 4: Analysis
% - calculate column integrated heat, salt and oxygen content
% - plot gradients in oxygen and show variability or "events"