function [x] = locateParticles(I, beadParameter)
% [x] = locateParticles(I, beadParameter) locates particles in the image
%
% INPUTS
% -------------------------------------------------------------------------
%   I:              Input volumetric image
%   beadParameter:  Parameters to detect particle position in images
%
% OUTPUTS
% -------------------------------------------------------------------------
%   x:              Voxel-level estimate of particle center in MxNxO format
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
thres = beadParameter.thres;    %Threshold value
minPixels = beadParameter.minSize;  %Minimum pixel count in blob for bead
maxPixels = beadParameter.maxSize;  %Maximum pixel count in blob for bead

% Image thresholding
BW = I>thres;

% Find bead blobs
CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
beadBlob = numPixels>minPixels & numPixels<maxPixels;

% Find centroid of all the eligible blob;
S = regionprops(CC,'Centroid');
blobPts = round(double(struct2dataset(S)));
blobPts = blobPts(beadBlob,:);
temp = blobPts;

% Convert to m,n,o coordinates
blobPts(:,1) = temp(:,2);
blobPts(:,2) = temp(:,1);
x = blobPts;

end

