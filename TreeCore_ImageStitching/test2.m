% clc
% close all
% clear all
% a= imread('Ar.jpg'); %read the left part of the image
% b=imread('Br.jpg');  %read the right part of the image
% sa= size(a); %get the size of the left image
% sb = size(b);%get the size of the left image
% b= imresize(b,[sa(1) sa(2)]); %now resize 'b' as per the size of 'a' in order to get perfect sized image
% c= [a b]; % club the image a & b into another new image after making their size same
% %note the missing semicolon in the above line inside bracket
% imshow(c) % show the image


original  = rgb2gray(imread('Br.jpg'));
imshow(original);
title('Base image');
% distorted = imresize(original,0.7); 
distorted = rgb2gray(imread('Ar.jpg'));
figure; imshow(distorted);
%title('Transformed image');

ptsOriginal  = detectSURFFeatures(original);
ptsDistorted = detectSURFFeatures(distorted);

[featuresOriginal,validPtsOriginal] = extractFeatures(original,ptsOriginal);
[featuresDistorted,validPtsDistorted] = extractFeatures(distorted,ptsDistorted);
index_pairs = matchFeatures(featuresOriginal,featuresDistorted);
matchedPtsOriginal  = validPtsOriginal(index_pairs(:,1));
matchedPtsDistorted = validPtsDistorted(index_pairs(:,2));
figure; 
showMatchedFeatures(original,distorted, matchedPtsOriginal,matchedPtsDistorted);7
%title('Matched SURF points,including outliers');
[tform,inlierPtsDistorted,inlierPtsOriginal] = estimateGeometricTransform(matchedPtsDistorted,matchedPtsOriginal, 'similarity');
figure; 
showMatchedFeatures(original,distorted, inlierPtsOriginal,inlierPtsDistorted);
%title('Matched inlier points');
outputView = imref2d(size(original));
Ir = imwarp(distorted,tform,'OutputView',outputView);
figure; imshow(Ir); 
title('Recovered image');

%https://www.mathworks.com/help/images/ref/imfuse.html
C = imfuse(original, Ir, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
figure; imshow(C); 
title('fuse image');

%Write result image to file.
imwrite(C, 'fused.png');