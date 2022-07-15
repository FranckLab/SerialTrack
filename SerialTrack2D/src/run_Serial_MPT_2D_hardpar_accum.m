% %%%%%%%%%%%%%%%%%% SerialTrack (2D cumulative mode) %%%%%%%%%%%%%%%%%
% Main file of code "SerialTrack"
% ***********************************************
% Dimension:            2D
% Particle rigidity:    hard 
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
% Date: 2020.12.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
%%%%% Load 2D images %%%%%
[file_name,Img,MPTPara_temp] = funReadImage2; close all;
disp('%%%%%% Load images: Done! %%%%%%'); fprintf('\n');

%%%%% Update MPTPara %%%%%
MPTPara.gridxyROIRange = MPTPara_temp.gridxyROIRange;
MPTPara.LoadImgMethod = MPTPara_temp.LoadImgMethod;
MPTPara.ImgSize = MPTPara_temp.ImgSize;

%%%%% Load image mask file(s) %%%%%
try
    if MaskFileLoadingMode == 1
        %%%%% Load one image mask file for all images %%%%%
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
%%%%% Particle detection parameters %%%%%
%%%%%%%%%% Atom dance %%%%%%%%%%%%%
%%%%% Bead Parameter %%%%%
% BeadPara.thres = 0.35 ;         % Threshold for detecting particles
% BeadPara.beadRad = 0;           % Estimated radius of a single particle [px]
% BeadPara.minSize = 2;           % Minimum area of a single particle [px^2]
% BeadPara.maxSize = 40;          % Maximum area of a single particle [px^2]
% BeadPara.winSize = [5, 5];      % Default (not used)
% BeadPara.dccd = [1,1];          % Default (not used)
% BeadPara.abc = [1,1];           % Default (not used)
% BeadPara.forloop = 1;           % Default (not used)
% BeadPara.randNoise = 1e-7;      % Default (not used)
% BeadPara.PSF = [];              % PSF function; Example: PSF = fspecial('disk', BeadPara.beadRad-1 ); % Disk blur
% BeadPara.color = 'black';       % Foreground (particle) color: options, 'white' or 'black'

%%%%%%%%%% Pipe %%%%%%%%%%%%%
%%%%% Bead Parameters %%%%%
% BeadPara.thres = 0.5;           % Threshold for detecting particles
% BeadPara.beadRad = 3;           % Estimated radius of a single particle [px]
% BeadPara.minSize = 2;           % Minimum area of a single particle [px^2]
% BeadPara.maxSize = 20;          % Maximum area of a single particle [px^2]
% BeadPara.winSize = [5, 5];      % Default (not used)
% BeadPara.dccd = [1,1];          % Default (not used)
% BeadPara.abc = [1,1];           % Default (not used)
% BeadPara.forloop = 1;           % Default (not used)
% BeadPara.randNoise = 1e-7;      % Default (not used)
% BeadPara.PSF = [];              % PSF function; Example: PSF = fspecial('disk', BeadPara.beadRad-1 ); % Disk blur
% BeadPara.color = 'white';       % Foreground (particle) color: options, 'white' or 'black'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ImgSeqNum = 1; % First reference image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Use the first uploaded mask file %%%%%
try
    if (MaskFileLoadingMode==1) || (MaskFileLoadingMode==3)
        currImg = Img{ImgSeqNum}.*MPTPara.ImgRefMask;
    elseif MaskFileLoadingMode==2
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

%%%%% If PSF is non-empty, SerialTrack performs deconvolution %%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Several methods to detect particles %%%%%
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
% % Not fully implemented yet
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
hold on; plot( parCoordA(:,1), parCoordA(:,2), 'bo');
view(2); box on; axis equal; axis tight; set(gca,'fontsize',18);
title('Detected particles in ref image','fontweight','normal');
 
%%%%% Report detected beads # %%%%%
disp(['Detected particle # in ref image: ',num2str(size(parCoordA,1))]);
disp('%%%%%% Detect particles: Done! %%%%%%'); fprintf('\n');
 


%% %%%%% Initialization %%%%%
%%%%%  MPT parameters %%%%%
% MPTPara.f_o_s = 30;              % Size of search field: max(|u|,|v|) [px]
% MPTPara.n_neighborsMax = 25;     % Max # of neighboring particles
% MPTPara.n_neighborsMin = 1;      % Min # of neighboring particles
% MPTPara.locSolver = 1;           % Local solver: 1-topology-based feature; 2-histogram-based feature first and then topology-based feature;
% MPTPara.gbSolver = 3;            % Global step solver: 1-moving least square fitting; 2-global regularization; 3-ADMM iterations
% MPTPara.smoothness = 1e-2;       % Coefficient of regularization
% MPTPara.outlrThres = 2;          % Threshold for removing outliers in TPT
% MPTPara.maxIterNum = 20;         % Max ADMM iteration number
% MPTPara.iterStopThres = 1e-2;    % ADMM iteration stopping threshold
% MPTPara.strain_n_neighbors = 20; % # of neighboring particles used in strain gauge
% MPTPara.strain_f_o_s = 60;       % Size of virtual strain gauge [px]
% MPTPara.usePrevResults = 0;      % Whether use previous results or not: 0-no; 1-yes;
% MPTPara.distMissing = 2;         % Distance threshold to check whether particle has a match or not [px]

%%%%%% To store results %%%%%
parCoord_prev = cell(length(Img),1);        parCoord_prev{1} = parCoordA;
track_A2B_prev = cell(length(Img)-1,1);     uv_B2A_prev = cell(length(Img)-1,1);
resultDisp = cell(length(Img)-1,1);         resultDefGrad = cell(length(Img)-1,1);
 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

track_ratio = nan*zeros(length(Img)-1,1);
DefType = 'exp'; defList = [2:1:length(Img)]';
fig_track_rat = figure; ax_TR = axes; hold on; plot(defList,track_ratio,'r^-.','linewidth',1);
% adjust_fig(fig_track_rat,ax_TR,'','',''); box on; 
title('Tracking ratio');
xlabel('Frame #'); ylabel('Incremental tracking ratio');
try axis([2,length(file_name),0,1]); catch, end
drawnow

for ImgSeqNum = 2 : length(Img)  % "ImgSeqNum" is the frame index
    
    disp(['====== Frame #',num2str(ImgSeqNum),' ======']);
    
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
    [parCoordB_temp,uv_B2A_temp,resultDisp_temp,resultDefGrad_temp,track_A2B_temp,~] = fun_SerialTrack_2D_HardPar( ...
        ImgSeqNum,defImg,BeadPara,MPTPara,parCoordA,parCoord_prev(2:end),uv_B2A_prev);
     
    %%%%% Store results %%%%%
    parCoord_prev{ImgSeqNum} = parCoordB_temp;
    uv_B2A_prev{ImgSeqNum-1} = uv_B2A_temp;  % cumulative displacement
    track_A2B_prev{ImgSeqNum-1} = track_A2B_temp;
    resultDisp{ImgSeqNum-1} = resultDisp_temp;
    resultDefGrad{ImgSeqNum-1} = resultDefGrad_temp;
          
    track_A2B = track_A2B_prev{ImgSeqNum-1}; 
    track_ratio(ImgSeqNum-1) = length(track_A2B(track_A2B>0))/size(parCoord_prev{1},1);
    
    plot(defList,track_ratio,'r^-.','linewidth',1);
    drawnow
    
end
  

%%%%% cumulative tracking ratio %%%%%
disp('%%%%% Calculate incremental tracking ratio %%%%%'); fprintf('\n');
track_ratio = zeros(length(Img)-1,1);
DefType = 'exp'; defList = [2:1:length(Img)]';
  
for ImgSeqNum = 2 : length(Img)
    track_A2B = track_A2B_prev{ImgSeqNum-1}; 
    track_ratio(ImgSeqNum-1) = length(track_A2B(track_A2B>0))/size(parCoord_prev{ImgSeqNum},1);      
end
 
fig=figure; ax=axes; hold on; plot(defList,track_ratio,'r^-.','linewidth',1);
set(gcf,'color','w'); box on; set(gca,'fontsize',18); 
% adjust_fig(fig,ax,'','','');
title('Tracking ratio','fontweight','normal');  
xlabel('Frame #'); % ylabel('Tracking ratio');
try axis([2,length(Img),0,1]); catch; end


%%%%% Save results %%%%%
disp('%%%%% SerialTrack hard particle tracking: Done! %%%%%'); fprintf('\n');
results_file_name = 'results_2D_hardpar.mat';
mkdir results
save(['./results/' results_file_name],'parCoord_prev','uv_B2A_prev','resultDisp','resultDefGrad','track_A2B_prev');
 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%% Visualize tracked cumulative displacement of each frame %%%%%
disp('%%%%% Plot tracked cumulative deformations %%%%%'); fprintf('\n');

%%%%% Experimental parameters %%%%%
try xstep = MPTPara.xstep; catch, xstep = 1; end % unit: um/px
try ystep = MPTPara.ystep; catch, ystep = xstep; end
try tstep = MPTPara.tstep; catch, tstep = 1; end % unit: us  
 
%%%%% Plot tracked incremental displacement field %%%%%
%%%%% Make a video %%%%%
v = VideoWriter('video_2D_accum.avi'); v.FrameRate = 5; open(v); figure,

for ImgSeqNum = 2:length(Img) % ImgSeqNum: Frame #
    
    % Displacement from tracked particles on deformed frame
    disp_A2B_parCoordB = -uv_B2A_prev{ImgSeqNum-1};
    parCoordB = parCoord_prev{ImgSeqNum};

    %%%%% Plot displacements %%%%%
    clf, plotCone2(parCoordB(:,1)*xstep,parCoordB(:,2)*ystep,disp_A2B_parCoordB(:,1)*xstep ,disp_A2B_parCoordB(:,2)*ystep );
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
if size(parCoord_prev,1) > 1  % Only for multiple frames

disp('%%%%% Merge tracked trajectory segments %%%%%'); fprintf('\n');

%%%%% Initialization %%%%%
resultDispCurr = resultDisp{1}; 
parCoordA = resultDispCurr.parCoordA; 
parCoordATraj = cell(size(parCoordA,1),1); 

%%%%% Compute and collect all trajectory segments %%%%%
for parInd = 1:size(parCoordA,1)
    
    for ImgSeqNum = 2:(size(resultDisp,1)+1)
    
        resultDispCurr = resultDisp{ImgSeqNum-1};
        parCoordB = resultDispCurr.parCoordB; 
        track_A2B = resultDispCurr.track_A2B; 
        
        if track_A2B(parInd) > 0
            parCoordATraj{parInd}(ImgSeqNum-1,1:2) = parCoordB(track_A2B(parInd),1:2);
        else
            parCoordATraj{parInd}(ImgSeqNum-1,1:2) = [nan,nan]; 
        end
    end
     
end
 

%%%%% Plot tracked trajectories %%%%% 
disp('%%%%% Plot tracked trajectories %%%%%'); fprintf('\n');
figure,
for parInd = 1:size(parCoordA,1)
    try
        wayPoints = parCoordATraj{parInd}; 
        if length(wayPoints(isnan(wayPoints(:,1))==1,1)) < 3
            if (size(resultDisp,1)+1)<4
                hold on; line(xstep*wayPoints(isnan(wayPoints(:,1))<1,1),ystep*wayPoints(isnan(wayPoints(:,1))<1,2),'linewidth',1.2); view(3); % straight lines
            else
                hold on; fnplt(cscvn(diag([xstep,ystep])*wayPoints(isnan(wayPoints(:,1))<1,:)'),'',1.);
            end
        end
        % if sum(1-isnan(wayPoints(:,1)))>1 % Don't show if there is only one point on the trajectory
        %     hold on; plot(xstep*wayPoints(:,1),ystep*wayPoints(:,2),'.','markersize',5);
        % end
    catch
    end
end

set(gca,'fontsize',18); view(2); box on; axis equal; axis tight;  
title('Tracked particle trajectory','fontweight','normal');
set(gca,'YDir','reverse'); xlabel('x'); ylabel('y');
axis(xstep*[MPTPara.gridxyROIRange.gridx(1), MPTPara.gridxyROIRange.gridx(2), ...
            MPTPara.gridxyROIRange.gridy(1), MPTPara.gridxyROIRange.gridy(2) ]);

end    % End of if:  Only work for multiple frames


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Press "Ctrl + C" and modify codes below to plot interpolated ' ...
      'displacements and strains on a uniform grid mesh']);
disp(['Press "Enter" key to keep running the code']);
pause; 


ImgSeqNum = 2; % TODO: assign a Frame #


%%%%% Previously tracked displacement field %%%%%
resultDispCurr = resultDisp{ImgSeqNum-1};
resultDefGradCurr = resultDefGrad{ImgSeqNum-1};
disp_A2B_parCoordB = resultDispCurr.disp_A2B_parCoordB;
parCoordB = resultDispCurr.parCoordB;
  
%%%%% Interpolate scatterred data to gridded data %%%%%
sxy = min([round(0.5*MPTPara.f_o_s),20])*[1,1]; % Step size for griddata
smoothness = 1e-3; % Smoothness for regularization; "smoothness=0" means no regularization

[x_Grid_refB,y_Grid_refB,u_Grid_refB]=funScatter2Grid2D(parCoordB(:,1),parCoordB(:,2),disp_A2B_parCoordB(:,1),sxy,smoothness);
[~,~,v_Grid_refB]=funScatter2Grid2D(parCoordB(:,1),parCoordB(:,2),disp_A2B_parCoordB(:,2),sxy,smoothness);

% Apply ROI image mask
[u_Grid_refB, v_Grid_refB] = funRmROIOutside(x_Grid_refB,y_Grid_refB,MPTPara.ImgRefMask,u_Grid_refB,v_Grid_refB);

% Build a displacement vector
uv_Grid_refB_Vector=[u_Grid_refB(:),v_Grid_refB(:)]'; uv_Grid_refB_Vector=uv_Grid_refB_Vector(:);

% Calculate deformation gradient
D_Grid = funDerivativeOp(size(x_Grid_refB,1),size(x_Grid_refB,2),mean(sxy)); % Central finite difference operator
F_Grid_refB_Vector=D_Grid*uv_Grid_refB_Vector; % {F}={D}{U}

% Change "uv_Grid_refB_Vector" and "F_Grid_refB_Vector" in physical world units
uv_Grid_refB_Vector_PhysWorld = [u_Grid_refB(:)*xstep,v_Grid_refB(:)*ystep]';
uv_Grid_refB_Vector_PhysWorld = uv_Grid_refB_Vector_PhysWorld(:);

F_Grid_refB_Vector_PhysWorld = [F_Grid_refB_Vector(1:4:end), ...             % F11
                                F_Grid_refB_Vector(2:4:end)*ystep/xstep, ... % F21
                                F_Grid_refB_Vector(3:4:end)*xstep/ystep, ... % F12
                                F_Grid_refB_Vector(4:4:end)]';               % F22
F_Grid_refB_Vector_PhysWorld = F_Grid_refB_Vector_PhysWorld(:);

%%%%% Cone plot grid data: displecement %%%%%
figure, plotCone2(x_Grid_refB*xstep,y_Grid_refB*ystep,u_Grid_refB*xstep,v_Grid_refB*ystep );
set(gca,'fontsize',18); view(2); box on; axis equal; axis tight; set(gca,'YDir','reverse');
title('Tracked displacement','fontweight','normal');
axis([xstep*MPTPara.gridxyROIRange.gridx(1), xstep*MPTPara.gridxyROIRange.gridx(2), ...
      ystep*MPTPara.gridxyROIRange.gridy(1), ystep*MPTPara.gridxyROIRange.gridy(2) ]);
  

%%%%% Generate an FE-mesh %%%%%
[coordinatesFEM_refB,elementsFEM_refB] = funMeshSetUp(x_Grid_refB*xstep,y_Grid_refB*ystep);

%%%%% Cone plot grid data: displacement %%%%%
Plotdisp_show(uv_Grid_refB_Vector_PhysWorld, coordinatesFEM_refB, elementsFEM_refB,[],'NoEdgeColor');
 
%%%%% Cone plot grid data: infinitesimal strain %%%%%
Plotstrain_show(F_Grid_refB_Vector_PhysWorld, coordinatesFEM_refB, elementsFEM_refB,[],'NoEdgeColor',xstep,tstep);
 

 










