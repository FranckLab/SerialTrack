function [x0] = locateBeads(I,beadParameter)

%Parameters
thres = beadParameter.thres;    %Threshold value
minPixels = beadParameter.minSize;  %Minimum pixel count in blob for bead
maxPixels = beadParameter.maxSize;  %Maximum pixel count in blob for bead

%Binary Image
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

% in m,n,o coordinates
blobPts(:,1) = temp(:,2);
blobPts(:,2) = temp(:,1);
x0 = blobPts;


end

