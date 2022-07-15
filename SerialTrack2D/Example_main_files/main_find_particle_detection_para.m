%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SerialTrack: to test chosen particle detection parameters
% ===================================================
% Author: Jin Yang, Ph.D.
% Email: jyang526@wisc.edu -or- aldicdvc@gmail.com 
% Date: 02/2022; 07/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
clear; close all; clc; clearvars -global
disp('*************************************************************');
disp('****** Welcome to SerialTrack particle detection test ******');
addpath( './function/','./src/');


%% 3.1.1 APP version
% run Trial_MPT_start.mlapp; 


%% 3.1.2 Script version

%%%%% Image binary mask file %%%%%
im_roi_mask_file_path = ''; % TODO: modify this line: e.g., '.\img_par2track_pipe\im_roi.mat';
 
%%%%% Load 2D images %%%%%
[file_name,Img,MPTPara_temp] = funReadImage2; close all;
disp('%%%%%% Load images: Done! %%%%%%'); fprintf('\n');

%%%%% Update MPTPara %%%%%
MPTPara.gridxyROIRange = MPTPara_temp.gridxyROIRange;
MPTPara.LoadImgMethod = MPTPara_temp.LoadImgMethod;
MPTPara.ImgSize = MPTPara_temp.ImgSize;

%%%%% Load image mask file %%%%%
try load(im_roi_mask_file_path); catch; end
try MPTPara.ImgRefMask = im_roi'; % Load stored image roi if existed
catch, MPTPara.ImgRefMask = ones(size(Img{1})); % Set up default image mask file
end
disp('%%%%%% Load image mask file: Done! %%%%%%'); fprintf('\n');


%% ====== Detect particles ======
%%%%% Particle detection parameters %%%%%
%%%%%%%%%% Pipe %%%%%%%%%%%%%
%%%%% Bead Parameter %%%%%
BeadPara.thres = 0.5;           % Threshold for detecting particles
BeadPara.beadRad = 3;           % Estimated radius of a single particle [px]
BeadPara.minSize = 2;           % Minimum area of a single particle [px^2]
BeadPara.maxSize = 20;          % Maximum area of a single particle [px^2]
BeadPara.winSize = [5, 5];      % By default
BeadPara.dccd = [1,1];          % By default
BeadPara.abc = [1,1];           % By default
BeadPara.forloop = 1;           % By default
BeadPara.randNoise = 1e-7;      % By default
BeadPara.PSF = [];              % PSF function; Example: PSF = fspecial('disk', BeadPara.beadRad-1 ); % Disk blur
BeadPara.color = 'white';       % Foreground (particle) color: options, 'white' or 'black'
if strcmp(BeadPara.color,'black')==1,  BeadPara.beadRad = 0; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ImgSeqNum = 1; % First reference image
currImg = Img{ImgSeqNum}.*MPTPara.ImgRefMask;
currImg = currImg(MPTPara.gridxyROIRange.gridx(1):MPTPara.gridxyROIRange.gridx(2), ...
                  MPTPara.gridxyROIRange.gridy(1):MPTPara.gridxyROIRange.gridy(2));

%%%%% If PSF is non-empty, perform deconvolution %%%%%
if ~isempty(BeadPara.PSF)
    currImg = deconvlucy(currImg,BeadPara.PSF,6);
    disp(['----- Deconvolution frame #',num2str(ImgSeqNum),' ------']);
end
%%%% figure, imshow(currImg);

 
%%%%% Pre-process bead image if bead color is black %%%%%
if strcmp(BeadPara.color,'black')
    ImgGauss = imgaussfilt(imgaussfilt(currImg,1),1); % figure, imshow(uint16(ImgGauss));
    ImgGauss(ImgGauss > BeadPara.thres*max(double(currImg(:)))) = 0;
    bw = imbinarize(uint16(ImgGauss),'adaptive','ForegroundPolarity','dark','Sensitivity',0.8); % figure, imshow(bws2);
    bws2 = bwareaopen(bw,round(pi*BeadPara.minSize^2)); % remove all object containing fewer than BeadPara.minSize
    removeobjradius = BeadPara.minSize; % fill a gap in the particles
    se = strel('disk',removeobjradius);
    bws2 = imclose(bws2,se);
    currImg2 = double(bws2); % figure, imshow(uint8(currImg2));
    currImg2_norm = double(currImg2)/max(double(currImg2(:)));
else
    currImg2_norm = double(currImg)/max(double(currImg(:)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Several methods to detect and localize particles %%%%%
%%%%% Method 1: TPT code %%%%%
% x{1}{ImgSeqNum} = locateBeads(double(currImg2)/max(double(currImg2(:))),BeadPara); % Detect particles
% x{1}{ImgSeqNum} = radial2center(double(currImg2)/max(double(currImg2(:))),x{1}{ImgSeqNum},BeadPara); % Localize particles
% ----------------------------
%%%%% Method 2: LoG operator (modified TracTrac code) %%%%%
%pre-process to get threshold and size values
YN = 0;
while YN == 0
    prompt = 'Adjust particle detection parameters? (Y/N) [N]: ';
    YN_ = input(prompt,'s');

    if strcmpi(YN_,'N')
        YN = 1;

    else

        BeadPara = funGetBeadPara(BeadPara,currImg2_norm);

        %run the particle detection and localization
        x{1}{ImgSeqNum} = f_detect_particles(currImg2_norm,BeadPara);
         
        %%%%% Store particle positions as "parCoordA" %%%%%
        x{1}{ImgSeqNum} = x{1}{ImgSeqNum} + [MPTPara.gridxyROIRange.gridx(1)-1, MPTPara.gridxyROIRange.gridy(1)-1];
        parCoordA = x{1}{ImgSeqNum};

        %%%%% Remove bad parCoord outside the image area %%%%%
        for tempi=1:2, parCoordA( parCoordA(:,tempi)>size(Img{ImgSeqNum},tempi), : ) = []; end
        for tempi=1:2, parCoordA( parCoordA(:,tempi)<1, : ) = []; end

        %%%%% Plot %%%%%
        figure, imshow(imread(file_name{1,1}))
        hold on; plot( parCoordA(:,1), parCoordA(:,2), 'r.');
        view(2); box on; axis equal; axis tight; set(gca,'fontsize',18);
        title('Detected particles in ref image','fontweight','normal');

    end

end
%%%%%%%%%%%% With decided particle detection parameters %%%%%%%%%%%
%run the particle detection and localization
x{1}{ImgSeqNum} = f_detect_particles(currImg2_norm,BeadPara);
% ----------------------------
%%%%% Method 3: sift code %%%%%
% [~,descriptors,locs] = sift(currImg2);
% x{1}{ImgSeqNum} = locs(:,1:2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Store particle positions as "parCoordA" %%%%%
x{1}{ImgSeqNum} = x{1}{ImgSeqNum} + [MPTPara.gridxyROIRange.gridx(1)-1, MPTPara.gridxyROIRange.gridy(1)-1];
parCoordA = x{1}{ImgSeqNum};

%%%%% Remove bad parCoord outside the image area %%%%%
for tempi=1:2, parCoordA( parCoordA(:,tempi)>size(Img{ImgSeqNum},tempi), : ) = []; end
for tempi=1:2, parCoordA( parCoordA(:,tempi)<1, : ) = []; end
 
%%%%% Plot %%%%%
figure, imshow(imread(file_name{1,1}))
hold on; plot( parCoordA(:,1), parCoordA(:,2), 'bo');
view(2); box on; axis equal; axis tight; set(gca,'fontsize',18);
title('Detected particles in ref image','fontweight','normal');
 
%%%%% Report detected beads # %%%%%
disp(['Detected particle # in ref image: ',num2str(size(parCoordA,1))]);
disp('%%%%%% Detect particles: Done! %%%%%%'); fprintf('\n');



