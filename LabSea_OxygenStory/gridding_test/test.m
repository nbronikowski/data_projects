clear; close all; clc; 

%%
% load data.mat % subset of data 
path_name = './../mat_files/'; var_name  = 'pearldiver_data';
load(fullfile(path_name,[var_name,'_oxy_qc.mat']));
dat = pearldiver; clear pearldiver;

%% 
idt = dat.dateNum>datenum(2020,02,01) &  dat.dateNum<datenum(2020,02,10);
X = dat.dateNum(idt);
Y = dat.depth(idt); 
Z = dat.adjusted_oxygen_concentration(idt);

Xi = X(1):datenum(seconds(60)):X(end);
Y = mean_interp(X,Y,Xi);
Z = mean_interp(X,Z,Xi);
X = Xi(:); clear Xi;
X = X-X(1);

% Get rid of nan's
idnan = isnan(Z);
X(idnan)=[];Y(idnan)=[];Z(idnan)=[];

% Define the grid    
[gridX, gridY] = meshgrid(0:0.5/24:max(X),0:5:1000);

epsilon = 0.5; k = 150;
% Z_grid = simpleGridData(X, Y, Z, gridX,gridY);
% Z_grid=RFBinterpolate(X,Y,Z,gridX,gridY,epsilon);
% Z_grid = knnIdwInterpolation(X, Y, Z, gridX, gridY, k);
% Z_grid = knnRbfInterpolation(X, Y, Z, gridX, gridY,1,100,0.0001);
% Z_grid = ordkriging(X, Y, Z, gridX, gridY);
% Z_grid = rbf(X, Y, Z, gridX, gridY, 'thin_plate_spline',1e-3);
% vfit = variogram([X(:), Y(:) ],Z(:),'plotit',true,'anisotropy',false);

Z_grid = barnes(X,Y,Z,gridX,gridY,1/4,20,5);

%% PLOT 
plot_result
