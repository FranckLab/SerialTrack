% %%%%%%%%%%%%%%%%%% SerialTrack (2D cumulative mode) %%%%%%%%%%%%%%%%%
% Main file of code "SerialTrack"
% ***********************************************
% Dimension:            2D
% Particle rigidity:    soft (particle shape is not rigid) 
% Tracking mode:        cumulative
% -----------------------------------------------
%
% -----------------------------------------------
% References
% [1] M Patel, SE Leggett, AK Landauer, IY Wong, C Franck. Rapid,
%     topology-based particle tracking for high-resolution measurements of
%     large complex 3D motion fields. Scientific Reports. 8:5581 (2018).
% [2] J Yang, L Hazlett, AK Landauer, C Franck. Augmented Lagrangian
%     Digital Volume Correlation (ALDVC). Experimental Mechanics (2020).
% [3] T Janke, R Schwarze, K Bauer. Part2Track: A MATLAB package for double
%     frame and time resolved Particle Tracking Velocimetry. 11, 100413, SoftwareX (2020).
% [4] J Heyman. TracTrac: a fast multi-object tracking algorithm for motion
%     estimation. Computers & Geosciences, vol 128, 11-18 (2019).
% [5] https://www.mathworks.com/matlabcentral/fileexchange/77347-gridded-interpolation-and-gradients-of-3d-scattered-data
% [6] https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend
% -----------------------------------------------
% Author: Jin Yang
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Date: 2020.12; 2022.07
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
%%%%% Load 2D images %%%%%
[file_name,Img,MPTPara_temp] = funReadImage2; close all;
disp('%%%%%% Load images: Done! %%%%%%'); fprintf('\n');

%%%%% Update MPTPara %%%%%
MPTPara.gridxyROIRange = MPTPara_temp.gridxyROIRange;
MPTPara.LoadImgMethod = MPTPara_temp.LoadImgMethod;
MPTPara.ImgSize = MPTPara_temp.ImgSize;

%%%%% Load image mask file %%%%%
try
    if MaskFileLoadingMode == 1
        %%%%% Load one image mask file for all %%%%%
        [file_mask_name,ImgMaskFiles] = funReadImageMask2; close all;
    elseif MaskFileLoadingMode == 2
        %%%%% Load multiple image mask files %%%%%
        [file_mask_name,ImgMaskFiles] = funReadImageMask2; close all;
    elseif MaskFileLoadingMode == 3
        try load(im_roi_mask_file_path); catch; end
        try MPTPara.ImgRefMask = im_roi'; % Load stored image roi if existed
        catch, MPTPara.ImgRefMask = ones(size(Img{1})); % Set up default image mask file
        end
    else
        MPTPara.ImgRefMask = ones(size(Img{1})); % Set up default image mask file
    end
    disp('%%%%%% Load image mask file: Done! %%%%%%'); fprintf('\n');
catch
    MPTPara.ImgRefMask = ones(size(Img{1})); % Set up default image mask file
end

   
%% ====== Detect and localize particles ======
%%%%% Particle detection and localization parameters %%%%%
% The user input "BeadPara"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ImgSeqNum = 1; % First reference image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Apply the uploaded mask file %%%%%
try
    if MaskFileLoadingMode == 1 || 3
        currImg = Img{ImgSeqNum}.*MPTPara.ImgRefMask;
    elseif MaskFileLoadingMode == 2
        currImg = Img{ImgSeqNum}.*ImgMaskFiles{ImgSeqNum};
    else
        currImg = Img{ImgSeqNum};  
    end
catch
    currImg = Img{ImgSeqNum};  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Apply the ROI region %%%%%
currImg = currImg(MPTPara.gridxyROIRange.gridx(1):MPTPara.gridxyROIRange.gridx(2), ...
                  MPTPara.gridxyROIRange.gridy(1):MPTPara.gridxyROIRange.gridy(2));
refImg = currImg; % Store first frame as "refImg"


%%%%% If PSF is non-empty, perform deconvolution %%%%%
if ~isempty(BeadPara.PSF)
    currImg = deconvlucy(currImg,BeadPara.PSF,6);
    disp(['----- Deconvolution frame #',num2str(ImgSeqNum),' ------']);
end
%%%% figure, imshow(currImg);

 
%%%%% Pre-process bead image if bead color is black %%%%%
if strcmp(BeadPara.color,'black')
%     ImgGauss = imgaussfilt(imgaussfilt(currImg,1),1); % figure, imshow(uint16(ImgGauss));
%     ImgGauss(ImgGauss > BeadPara.thres*max(double(currImg(:)))) = 0;
%     bw = imbinarize(uint16(ImgGauss),'adaptive','ForegroundPolarity','dark','Sensitivity',0.8); % figure, imshow(bws2);
%     bws2 = bwareaopen(bw,BeadPara.minSize); % remove all object containing fewer than BeadPara.minSize
%     removeobjradius = sqrt(BeadPara.minSize/pi); % fill a gaps in particles
%     se = strel('disk',round(removeobjradius));
%     bws2 = imclose(bws2,se);
    currImg_norm = double(currImg)/max(double(currImg(:)));
    currImg2_norm = imcomplement(currImg_norm); % figure, imshow(uint8(currImg2));
else
    currImg2_norm = double(currImg)/max(double(currImg(:)));
end

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
% [~,descriptors,locs] = sift(currImg);
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
hold on; plot( parCoordA(:,1), parCoordA(:,2), 'r.');
view(2); box on; axis equal; axis tight; set(gca,'fontsize',18);
title('Detected particles in ref image','fontweight','normal');
 
%%%%% Report detected beads # %%%%%
disp(['Detected particle # in ref image: ',num2str(size(parCoordA,1))]);
disp('%%%%%% Detect particles: Done! %%%%%%'); fprintf('\n');
 

 
%% %%%%% Initialization %%%%%
%%%%%  MPT Parameter %%%%%
% The user input MPTPara

%%%%%% To store results %%%%%
xyGrid_prev = cell(length(Img)-1,1);
uvGrid_B2A_prev = cell(length(Img)-1,1);
track_A2B_prev = cell(length(Img)-1,1);
parCoord_prev = cell(length(Img),1); parCoord_prev{1} = parCoordA;

 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
fig=figure; ax=axes; defList = [2:1:length(Img)]'; track_ratio = nan*defList;
hold on; plot(defList,track_ratio,'r^-.','linewidth',1);
adjust_fig(fig,ax,'','',''); box on; title('');
xlabel('Frame #'); ylabel('Tracking ratio');
axis([2,length(Img),0,1]);
drawnow

for ImgSeqNum = 2 : length(Img)  % "ImgSeqNum" is the frame index
    
    close all; disp(['====== Frame #',num2str(ImgSeqNum),' ======']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Apply the uploaded mask file %%%%%
    clear defImg;
    try 
        if (MaskFileLoadingMode==1) || (MaskFileLoadingMode==3)
            defImg = Img{ImgSeqNum}.*MPTPara.ImgRefMask;
        elseif MaskFileLoadingMode==2
            defImg = Img{1,ImgSeqNum}.*ImgMaskFiles{1,ImgSeqNum};
        else
            defImg = Img{ImgSeqNum};
        end
    catch
       defImg = Img{ImgSeqNum};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% SerialTrack particle tracking %%%%%
    [temp1,temp2,track_A2B,parCoordB] = fun_SerialTrack_2D_SoftPar( ...
       ImgSeqNum,defImg,BeadPara,MPTPara,parCoord_prev{1},xyGrid_prev,uvGrid_B2A_prev);
      
    %%%%% Store results %%%%%
    xyGrid_prev{ImgSeqNum-1} = temp1;        
    uvGrid_B2A_prev{ImgSeqNum-1} = temp2;
    track_A2B_prev{ImgSeqNum-1} = track_A2B;
    parCoord_prev{ImgSeqNum} = parCoordB;
    
    track_A2B = track_A2B_prev{ImgSeqNum-1}; 
    track_ratio(ImgSeqNum-1) = length(track_A2B(track_A2B>0))/size(parCoord_prev{1},1);
    
    plot(defList,track_ratio,'r^-.','linewidth',1);
    drawnow
      
   
end
  
%% %%%%% cumulative tracking ratio %%%%%
disp('%%%%% Calculate cumulative tracking ratio %%%%%'); fprintf('\n');
track_ratio = zeros(length(Img)-1,1);
DefType = 'exp'; defList = [2:1:length(Img)]';
  
for ImgSeqNum = 2 : length(Img)
    track_A2B = track_A2B_prev{ImgSeqNum-1}; 
    track_ratio(ImgSeqNum-1) = length(track_A2B(track_A2B>0))/size(parCoord_prev{1},1);     
end
 
fig=figure; ax=axes; hold on; plot(defList,track_ratio,'r^-.','linewidth',1);
adjust_fig(fig,ax,'','',''); box on; title('');
xlabel('Frame #'); ylabel('Tracking ratio');
axis([2,length(Img),0,1]);


%%%%% Save results %%%%%
disp('===== SerialTrack soft particle tracking: Done! =====');  
results_file_name = 'results_2D_softpar.mat';
mkdir results
save(['./results/' results_file_name],'xyGrid_prev','uvGrid_B2A_prev','track_A2B_prev','parCoord_prev');
 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%% Visualize tracked cumulative displacement of each frame %%%%%
disp('%%%%% Plot tracked cumulative deformations %%%%%'); fprintf('\n');

%%%%% Experimental parameters %%%%%
try xstep = MPTPara.xstep; catch, xstep = 1; end % unit: um/px
try tstep = MPTPara.tstep; catch, tstep = 1; end % unit: us  
% ImgSeqNum  % Frame #
 
%%%%% Plot tracked incremental displacement field %%%%%
%%%%% Make a video %%%%%
v = VideoWriter('video_2D_accum_softpar.avi'); v.FrameRate = 5; open(v); figure,

for ImgSeqNum = [ 51 ] % 2 : length(Img) % ImgSeqNum: Frame #
    
    % Displacement from tracked particles on deformed frame
    disp_A2B = -uvGrid_B2A_prev{ImgSeqNum-1};
    xy_refA = xyGrid_prev{ImgSeqNum-1};
    parCoordB = parCoord_prev{ImgSeqNum};
    track_A2B = track_A2B_prev{ImgSeqNum-1};

    %%%%% Interpolate to scattered points %%%%%
    Fx = scatteredInterpolant(xy_refA,disp_A2B(:,1),'linear','linear');
    disp_A2B_uScatter = Fx(parCoordA(track_A2B>0,1:2));
    Fy = scatteredInterpolant(xy_refA,disp_A2B(:,2),'linear','linear');
    disp_A2B_vScatter = Fy(parCoordA(track_A2B>0,1:2));
    
    xScatter_refA = parCoordA((track_A2B>0),1);
    yScatter_refA = parCoordA((track_A2B>0),2);
    xScatter_refB = parCoordB(track_A2B(track_A2B>0),1);
    yScatter_refB = parCoordB(track_A2B(track_A2B>0),2);
    
    % %%%%% Plot displacements refA %%%%%
    % clf,  plotCone2(xScatter_refA*xstep,yScatter_refA*xstep,disp_A2B_uScatter*xstep,disp_A2B_vScatter*xstep);
    % set(gca,'fontsize',18); view(2); box on; axis equal; axis tight; set(gca,'YDir','reverse');
    % title(['Tracked cumulative displacement (#',num2str(ImgSeqNum),')'],'fontweight','normal');
    % xlabel(''); ylabel('');
    
    %%%%% Plot displacements refB %%%%%
    clf,  plotCone2(xScatter_refB*xstep,yScatter_refB*xstep,disp_A2B_uScatter*xstep,disp_A2B_vScatter*xstep);
    set(gca,'fontsize',18); view(2); box on; axis equal; axis tight; set(gca,'YDir','reverse');
    title(['Tracked cumulative disp (#',num2str(ImgSeqNum),')'],'fontweight','normal');
    xlabel('x'); ylabel('y');
    axis(xstep*[MPTPara.gridxyROIRange.gridx(1), MPTPara.gridxyROIRange.gridx(2), ...
                MPTPara.gridxyROIRange.gridy(1), MPTPara.gridxyROIRange.gridy(2) ]);
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    
end
close(v);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%% Compute trajectory %%%%%

%%%%% Initialization %%%%%
parCoordATraj = cell(size(parCoordA,1),1); 

%%%%% Compute and collect all trajectory segments %%%%%
for parInd = 1:size(parCoordA,1)
    
    for ImgSeqNum = 2 : (size(parCoord_prev,1))
    
        parCoordB = parCoord_prev{ImgSeqNum};
        track_A2B = track_A2B_prev{ImgSeqNum-1}; 
        try
        if track_A2B(parInd) > 0
            parCoordATraj{parInd}(ImgSeqNum-1,1:2) = parCoordB(track_A2B(parInd),1:2);
        else
            parCoordATraj{parInd}(ImgSeqNum-1,1:2) = [nan,nan]; 
        end
        catch
        end
    end
     
end
 

%%%%% Plot tracked trajectories %%%%% 
disp('%%%%% Plot tracked trajectories %%%%%'); fprintf('\n');
figure,
for parInd = 1:size(parCoordA,1)
    try
        wayPoints = parCoordATraj{parInd}; 
        if length(wayPoints(isnan(wayPoints(:,1))==1,1)) < size(parCoordA,1)-minTrajSegLength 
            if (size(parCoord_prev,1))<4
                hold on; line(xstep*wayPoints(isnan(wayPoints(:,1))<1,1),xstep*wayPoints(isnan(wayPoints(:,1))<1,2),'linewidth',1.2);  % straight lines
            else
                hold on; fnplt(cscvn(xstep*wayPoints(isnan(wayPoints(:,1))<1,:)'),'',1.2);
            end
        end
        % if sum(1-isnan(wayPoints(:,1)))>1 % Don't show if there is only one point on the trajectory
        %     hold on; plot(xstep*wayPoints(:,1),xstep*wayPoints(:,2),'.','markersize',5);
        % end
    catch
    end
end

set(gca,'fontsize',18); view(2); box on; axis equal; axis tight;  
title('Tracked particle trajectory','fontweight','normal');
set(gca,'YDir','reverse'); xlabel(''); ylabel('');
axis(xstep*[MPTPara.gridxyROIRange.gridx(1), MPTPara.gridxyROIRange.gridx(2), ...
            MPTPara.gridxyROIRange.gridy(1), MPTPara.gridxyROIRange.gridy(2) ]);


        
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Press "Ctrl + C" and modify codes below to plot interpolated displacements and strains on a uniform grid mesh');
disp(['Press "Enter" key to keep running the code']);
pause;

ImgSeqNum = 4; % Frame #


%%%%% Previously tracked displacement field %%%%%
disp_A2B = -uvGrid_B2A_prev{ImgSeqNum-1};
xy_refA = xyGrid_prev{ImgSeqNum-1};
parCoordB = parCoord_prev{ImgSeqNum}; 

sxy = min([round(0.5*MPTPara.f_o_s),20])*[1,1];  % Please don't change sxy
[xGrid_refA,yGrid_refA] = ndgrid(MPTPara.gridxyROIRange.gridx(1) : sxy(1) : MPTPara.gridxyROIRange.gridx(2), ...
                                 MPTPara.gridxyROIRange.gridy(1) : sxy(2) : MPTPara.gridxyROIRange.gridy(2));
      
%%%%% TODO: JY!!!
% Apply ROI image mask
% [u_Grid_refB, v_Grid_refB] = funRmROIOutside(x_Grid_refB,y_Grid_refB,MPTPara.ImgRefMask,u_Grid_refB,v_Grid_refB);

%%%%% Compute F, or deformation gradient tensor on deformed frame %%%%%
[M,N] = size(xGrid_refA);
DMat = funDerivativeOp(M,N,mean(sxy)); % Central finite difference operator
uVec = [disp_A2B(:,1), disp_A2B(:,2)]'; uVec=uVec(:);
FVec = DMat*uVec; % {F}={D}{U}

%%%%% Generate an FE-mesh %%%%%
[coordinatesFEM_refA,elementsFEM_refA] = funMeshSetUp(xGrid_refA*xstep,yGrid_refA*xstep);

%%%%% Cone plot grid data: displacement %%%%%
Plotdisp_show(uVec*xstep, coordinatesFEM_refA*xstep, elementsFEM_refA, [],'NoEdgeColor');
 
%%%%% Cone plot grid data: infinitesimal strain %%%%%
Plotstrain_show(FVec,coordinatesFEM_refA*xstep,elementsFEM_refA, [],'NoEdgeColor');

 










