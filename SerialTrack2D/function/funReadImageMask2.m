function [file_mask_name,ImgMaskFile] = funReadImageMask2(varargin)
%FUNCTION [file_name,Img,TPTpara] = ReadImage2(varargin)
% MATLAB script: ReadImage.m
% ----------------------------------------------
%   This script is to load DIC images 
%   Images can be loaded by:
%       i) selecting a folder which included all the DIC raw images, 
%       ii) inputing image file name prefix keywords
%       iii) manually select DIC raw images
%
%   INPUT: No inputs are needed
%
%   OUTPUT: file_name    Loaded DIC raw image file name
%           Img          Loaded DIC images
%           TPTpara      DIC parameters
%
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 02/2020.
% ==============================================

%%
fprintf('Choose method to load image mask files:  \n')
fprintf('     0: Select image mask file folder;  \n')
fprintf('     1: Use prefix of image mask file names;  \n')
fprintf('     2: Manually select image mask files.  \n')
prompt = 'Input here: ';
LoadImgMethod = input(prompt);

switch LoadImgMethod 
    case 0
        % ==============================================
        imgfoldername = uigetdir(pwd,'Select image mask file folder');
        addpath([imgfoldername,'\']);
        img1 = dir(fullfile(imgfoldername,'*.jpg'));
        img2 = dir(fullfile(imgfoldername,'*.jpeg'));
        img3 = dir(fullfile(imgfoldername,'*.tif'));
        img4 = dir(fullfile(imgfoldername,'*.tiff'));
        img5 = dir(fullfile(imgfoldername,'*.bmp'));
        img6 = dir(fullfile(imgfoldername,'*.png'));
        img7 = dir(fullfile(imgfoldername,'*.jp2'));
        file_mask_name = [img1;img2;img3;img4;img5;img6;img7];
        file_mask_name = struct2cell(file_mask_name);
    case 1
        % ==============================================
        fprintf('What is prefix of DIC image mask files? E.g. img_0*.tif.   \n')
        prompt = 'Input here: ';
        file_mask_name = input(prompt,'s');
        [~,imgname,imgext] = fileparts(file_mask_name);
        file_mask_name = dir([imgname,imgext]);
        file_mask_name = struct2cell(file_mask_name);
    otherwise
        % ==============================================
        disp('--- Please load first image mask file ---')
        file_mask_name{1,1} = uigetfile('*.tif','Select reference image mask file');
        % disp('--- Please load next image mask file ---')
        % file_mask_name{1,2} = uigetfile('*.tif','Select deformed image mask file');
        prompt = 'Do you want to load more deformed image mask files? (0-Yes; 1-No)';
        DoYouWantToLoadMoreImages = input(prompt); imageNo = 1;
        while ( DoYouWantToLoadMoreImages == 0 )   
            imageNo = imageNo + 1;
            file_mask_name{1,imageNo} = uigetfile('*.tif','Select deformed image mask file');
            prompt = 'Do you want to load more deformed image mask files? (0-Yes; 1-No)';
            DoYouWantToLoadMoreImages = input(prompt);
        end
end

% ==============================================
% The following codes only consider two images comparasion
numImages = size(file_mask_name,2);
for i = 1:numImages
    ImgMaskFile{i} = imread(file_mask_name{1,i});
    % Change color RGB images to grayscale images
    [~, ~, numberOfColorChannels] = size(ImgMaskFile{i});
    if (numberOfColorChannels==3)
        ImgMaskFile{i} = rgb2gray(ImgMaskFile{i});
    end
    temp1 = double((ImgMaskFile{i}))';
    temp1(~temp1) = nan;
    ImgMaskFile{i} = temp1;
end

    

end

