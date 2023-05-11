
clear; clc; close all
S=readgeotable('job_LAM12092019_corrected_14092019_MERGE.shp');
 
X = S.Crrctd_X;
Y = S.Crrctd_Y;
Z = S.Z_Elevatio;
pp = double(S.Point);

idx = 333;
Za = Z;
% Za(333:end) = flipud(Za(333:end));
Za(333:end) = Za(333:end)-0.2;
TF = isoutlier(Za) ;
% id1 = find(X==1.023083000000000e+03 & Y==9.324720000000000e+03);
id11= find(X==1.016840383320000e+03 & Y==9.273267338000000e+02);
% id2 = find(X==1.016840383320000e+03 & Y==9.273267338000000e+03);
id22= find(X==1.023083000000000e+03 & Y==9.324720000000000e+02);
id3 = find(X==1.015938837640000e+03 & Y==9.293073792480000e+02);
id4 = find(X==1.023988234840000e+03 & Y==9.290105801640000e+02);
% Za(id1)=NaN;
% Za(id2)=NaN;
Za(id3:id4)=NaN;
Za(id11)=NaN;
Za(id22)=NaN;

[tmpX1,~,tmpX3] = unique(X);
[tmpY1,~,tmpY3] = unique(Y);
[Xi,Yi] = meshgrid(tmpX1,tmpY1);
tmpZ = sub2ind(size(Xi),tmpY3,tmpX3); % use the third output from unique as indices
Zi = NaN(size(Xi)); % preallocate Z
Zi(tmpZ) = Za; 

Zin =inpaint_nans(Zi);

figure()
subplot(121)
scatter(X,Y,100,Z,'filled');
colormap(cmocean('amp'))

subplot(122)
scatter(X,Y,100,Za,'filled');
colormap(cmocean('amp'))



figure()
imagescn(Xi,Yi,smooth2a(Zin,2,2)); hold on
scatter(X,Y,50,Za,'filled');
colorbar;
colormap(cmocean('amp'))

figure()
surf(Xi,Yi,smooth2a(Zin,2,2))
shading flat;  colorbar; colormap(cmocean('amp'))
title('Pauls Depression')






% myMeanFunction = @(block_struct) nanmean(nanmean(block_struct.data));
% ZNaNMeans = blockproc(Zi, [2, 2], myMeanFunction);
% XNaNMeans = blockproc(Xi, [2, 2], myMeanFunction);
% YNaNMeans = blockproc(Yi, [2, 2], myMeanFunction);

% 
% 
% 
% r1x = nanmean(diff(X));  % small influence
% r2x = nanmean(diff(X))*150;  % large cutoff
% r1y = nanmean(diff(X));   % same but for y
% r2y = nanmean(diff(Y))*50;   % 
% mtd = 1;     % exponential method, 2 = gaussian
% snr = 0.666; % signal to noise ratio
% N  = 200;    % max size of local cov. matrix
% N_r=  20;    % num. rand. datapoints
% % 
% 
% Zn = true_obana(Z,X,Y,Xi,Yi,...
%     r1x,r2x,r1y,r2y,mtd,snr,N,N_r);
% 
% figure(); hold on
% 
% imagescn(Xi,Yi,Zn)
% scatter(X,Y,50,Z   ,'filled')
% 
% colormap(cmocean('amp'))
% 
% colorbar
