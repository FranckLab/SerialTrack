function BeadPara = funGetBeadPara(BeadPara,Img)
% MATLAB script: funGetBeadPara.m
% ----------------------------------------------
%   This script is used to collect Bead Parameter settings from the user
%   Careful selection of parameters is important to correctly segment
%   particles from the background, by specifying the threshold and size
%   parameters
%
%   INPUT: BeadPara      Struct of bead parameters to use for segmentation
%          Img           An example input image (usually the first)
%
%   OUTPUT: BeadPara    Updated bead parameters struct
%
% ----------------------------------------------
% Author: Alex Landauer
% Last time updated: 07/2022.
% ==============================================


% first, get the threshold estimate
f1 = figure;
histogram(Img(Img~=0))
drawnow
f2 = figure;
imshow(Img,[]),colorbar,axis image
drawnow

BeadPara.thres = input('Enter binarization threshold estimate: ');
try
    close(f1)
    close(f2)
catch
end
YN = 0;
while YN == 0
    f3 = figure;
    imshow(Img>BeadPara.thres,[]),colorbar,axis image
    drawnow
    
    YN_ = input('Binarization okay? (Y/N) [N]: ','s');
    
    if strcmpi(YN_,'Y')
        YN = 1;
    else
        f1 = figure;
        histogram(Img(Img~=0))
        drawnow
        f2 = figure;
        imshow(Img,[]),colorbar,axis image
        drawnow
        BeadPara.thres = input('Enter binarization threshold estimate from histogram: ');
    end
    
end
try
    close(f1)
catch
end
try
    close(f2)
catch
end
try
    close(f3)
catch
end


% second, using this threshhold get the min and max sizes
% Parameters
thres = BeadPara.thres;    %Threshold value

% other params, set by user set later in this method
% minPixels = beadParameter.minSize;  %Minimum pixel count in blob for bead
% maxPixels = beadParameter.maxSize;  %Maximum pixel count in blob for bead

disp('%%%%%% Starting Binarization %%%%%%')

%Binary Image % figure, imshow(im);
BW = Img>thres; % figure, imshow(BW);

% Find bead blobs
CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);

%check the histogram of sizes and get user defined size limits
try
    nbins = sshist(numPixels);
catch
    nbins = 1;
end
nbins = max([15,nbins]);

h = figure;
histogram(numPixels,nbins)
disp('Set estimated bead rad as "0" if bead radii are various.')
estimatedRad = input('Enter estimated bead rad (px): ');  %Estimated blob radius for bead
minPixels = input('Enter min bead size (area): ');  %Minimum pixel count in blob for bead
maxPixels = input('Enter max bead size (area): ');  %Maximum pixel count in blob for bead
try
    close(h)
catch
end

BeadPara.beadRad = estimatedRad; %Estimated blob radius for bead
BeadPara.minSize = minPixels;  %Minimum pixel count in blob for bead
BeadPara.maxSize = maxPixels;  %Maximum pixel count in blob for bead




