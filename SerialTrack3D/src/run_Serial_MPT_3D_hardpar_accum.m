% %%%%%%%%%%%%%%%%%% SerialTrack (3D cumulative mode) %%%%%%%%%%%%%%%%%
% Main file of code "SerialTrack"
% ***********************************************
% Dimension:            3D
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
%%%%% Load 3D volumetric images %%%%%
try if isempty(fileFolder)~=1, cd(fileFolder); end; catch; end % Open image folder

ImgSeqNum=1; [file_name,Img] = funReadImage3(fileNameAll,ImgSeqNum); % Load image

try if isempty(fileFolder)~=1, cd(SerialTrackPath); end; catch; end % Come back to the main path

%%%%% Update MPTPara %%%%%
MPTPara.gridxyzROIRange.gridx = [1,size(Img{1},1)];
MPTPara.gridxyzROIRange.gridy = [1,size(Img{1},2)];
MPTPara.gridxyzROIRange.gridz = [1,size(Img{1},3)];

% figure, imagesc3(Img{1}) % To see volumetric image
disp('%%%%%% Load reference image: Done! %%%%%%'); fprintf('\n');

%%%%% Load image mask file %%%%%
try load(im_roi_mask_file_path); catch; end
try MPTPara.ImgRefMask = im_roi'; % Load stored image roi if existed
catch, MPTPara.ImgRefMask = ones(size(Img{1})); % Set up default image mask file
end
disp('%%%%%% Load image mask file: Done! %%%%%%'); fprintf('\n');

   
%% ====== Detect particles ======
%%%%% Particle detection parameters %%%%%
%%%%% Bead Parameters %%%%%
% BeadPara.detectionMethod = 2;   % Particle detection method: 1 = TPT (blob finding + radial projection), 
%                                                              2 = TracTrac (LoG blob finding + lsq fit of gaussian)
% BeadPara.thres = 0.4;           % Threshold for detecting particles
% BeadPara.beadRad = 0;           % Estimated radius of a single particle [px]
% BeadPara.minSize = 2;           % Minimum volume of a single particle [px^3]
% BeadPara.maxSize = 1000;        % Maximum volume of a single particle [px^3]
% BeadPara.winSize = [5, 5, 5];   % Default window search size for a single particle [px]
% BeadPara.dccd = [1,1,1];        % Default [not currently used]
% BeadPara.abc = [1,1,1];         % Default [not currently used]
% BeadPara.forloop = 1;           % Default [not currently used]
% BeadPara.randNoise = 1e-7;      % Default [not currently used]
% BeadPara.PSF = [];              % PSF function; Example: PSF = fspecial('disk', BeadPara.beadRad-1 ); % Disk blur
% BeadPara.color = 'white';       % Foreground (particle) color: options, 'white' or 'black'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ImgSeqNum = 1; % First reference image
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Several methods to detect particles %%%%%
try 
    BeadPara.detectionMethod = BeadPara.detectionMethod;
catch
    BeadPara.detectionMethod = 2;
end
%%%%% Method 1: TPT code %%%%%
if BeadPara.detectionMethod == 1 
    x{1}{ImgSeqNum} = locateParticles(double(Img{ImgSeqNum})/max(double(Img{ImgSeqNum}(:))),BeadPara); % Detect particles
    x{1}{ImgSeqNum} = radialcenter3dvec(double(Img{ImgSeqNum}),x{1}{ImgSeqNum},BeadPara); % Localize particles
% ----------------------------
%%%%% Method 2: Modified TracTrac code %%%%%
elseif BeadPara.detectionMethod == 2
    x{1}{ImgSeqNum} = f_detect_particles3(double(Img{ImgSeqNum})/max(double(Img{ImgSeqNum}(:))),BeadPara);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Store particle positions as "parCoordA" %%%%%
x{1}{ImgSeqNum} = x{1}{ImgSeqNum} + [MPTPara.gridxyzROIRange.gridx(1)-1, MPTPara.gridxyzROIRange.gridy(1)-1, MPTPara.gridxyzROIRange.gridz(1)-1];
parCoordA = x{1}{ImgSeqNum};

%%%%% Remove bad parCoord outside the image area %%%%%
for tempi=1:3, parCoordA( parCoordA(:,tempi)>size(Img{ImgSeqNum},tempi), : ) = []; end
for tempi=1:3, parCoordA( parCoordA(:,tempi)<1, : ) = []; end
 
%%%%% Plot %%%%%
% figure, plot3(parCoordA(:,1),parCoordA(:,2),parCoordA(:,3),'bo');
% view(3); box on; axis equal; axis tight; set(gca,'fontsize',18); 
% title('Detected particles in ref image','fontweight','normal');
 
%%%%% Report detected beads # %%%%%
disp(['Detected particle # in ref image: ',num2str(size(parCoordA,1))]);
disp('%%%%%% Detect particles: Done! %%%%%%'); fprintf('\n');
 


%% %%%%% Initialization %%%%%
%%%%% MPT Parameters %%%%%
% MPTPara.f_o_s = 60;              % Size of search field: max(|u|,|v|,|w|) [px]
% MPTPara.n_neighborsMax = 25;     % Max # of neighboring particles
% MPTPara.n_neighborsMin = 1;      % Min # of neighboring particles
% MPTPara.gbSolver = 2;            % Global step solver: 1-moving least square fitting; 2-global regularization; 3-ADMM iterations
% MPTPara.smoothness = 1e-1;       % Coefficient of regularization
% MPTPara.outlrThres = 5;          % Threshold for removing outliers in TPT
% MPTPara.maxIterNum = 20;         % Max ADMM iteration number
% MPTPara.iterStopThres = 1e-2;    % ADMM iteration stopping threshold
% MPTPara.strain_n_neighbors = 20; % # of neighboring particles used in strain gauge
% MPTPara.strain_f_o_s = 60;       % Size of virtual strain gauge [px]
% MPTPara.usePrevResults = 0;      % Whether use previous results or not: 0-no; 1-yes;
% MPTPara.distMissing = 5;         % Distance threshold to check whether particle has a match or not [px]


%%%%%% To store results %%%%%
parCoord_prev = cell(length(file_name),1);     parCoord_prev{1} = parCoordA;
uvw_B2A_prev = cell(length(file_name)-1,1);    track_A2B_prev = cell(length(file_name)-1,1);
resultDisp = cell(length(file_name)-1,1);      resultDefGrad = cell(length(file_name)-1,1);
 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ImgSeqNum = 2 : length(file_name)  % "ImgSeqNum" is the frame index
    
    disp(['====== Frame #',num2str(ImgSeqNum),' ======']);

    %%%%% Load image volumetric data %%%%%
    try if isempty(fileFolder)~=1, cd(fileFolder); end; catch; end % Open image folder
    tempvol = load(file_name{ImgSeqNum}); fieldName = fieldnames(tempvol);
    Img{2} = getfield(tempvol,fieldName{1}); clear tempvol; %#ok<GFLD>
    if iscell(Img{2}), Img{2}=Img{2}{1}; end
    try if isempty(fileFolder)~=1, cd(SerialTrackPath); end; catch; end % Come back to the main path

    %%%%% SerialTrack particle tracking %%%%%
    [parCoordB_temp,uvw_B2A_temp,resultDisp_temp,resultDefGrad_temp,track_A2B_temp,~] = fun_SerialTrack_3D_HardPar( ...
       ImgSeqNum,Img{2},BeadPara,MPTPara,parCoordA,parCoord_prev(2:end),uvw_B2A_prev);
    
    %%%%% Store results %%%%%
    parCoord_prev{ImgSeqNum} = parCoordB_temp;
    uvw_B2A_prev{ImgSeqNum-1} = uvw_B2A_temp;  % cumulative displacement
    resultDisp{ImgSeqNum-1} = resultDisp_temp;
    resultDefGrad{ImgSeqNum-1} = resultDefGrad_temp;
    track_A2B_prev{ImgSeqNum-1} = track_A2B_temp;
      
end
  

%% %%%%% cumulative tracking ratio %%%%%
disp('%%%%% Calculate cumulative tracking ratio %%%%%'); fprintf('\n');
track_ratio = zeros(length(file_name)-1,1);
DefType = 'exp'; defList = [2:1:length(file_name)]';
  
for ImgSeqNum = 2 : length(file_name)
    track_A2B = track_A2B_prev{ImgSeqNum-1}; 
    track_ratio(ImgSeqNum-1) = length(track_A2B(track_A2B>0))/size(parCoord_prev{ImgSeqNum},1);      
end
 
fig=figure; ax=axes; hold on; plot(defList,track_ratio,'r^-.','linewidth',1);
adjust_fig(fig,ax,'','',''); box on; title('');
xlabel('Frame #'); ylabel('Tracking ratio');
try axis([2,length(file_name),0,1]); catch; end

%%%%% Save results %%%%%
disp('%%%%%% SerialTrack 3D hard particle tracking: Done! %%%%%%'); fprintf('\n');
results_file_name = 'results_3D_hardpar.mat';
mkdir results
save(['./results/' results_file_name],'parCoord_prev','uvw_B2A_prev','resultDisp','resultDefGrad','track_A2B_prev');
 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%% Visualize tracked cumulative displacement of each frame %%%%%
disp('%%%%% Plot tracked cumulative deformations %%%%%'); fprintf('\n');

%%%%% Experimental parameters %%%%%
try xstep = MPTPara.xstep; catch, xstep = 1; end % unit: um/px
try ystep = MPTPara.ystep; catch, ystep = xstep; end
try zstep = MPTPara.zstep; catch, zstep = xstep; end
try tstep = MPTPara.tstep; catch, tstep = 1; end % unit: us  
 
%%%%% Plot tracked incremental displacement field %%%%%
%%%%% Make a video %%%%%
v = VideoWriter('video_3D_accum.avi'); v.FrameRate = 5; open(v); figure,

for ImgSeqNum = 2:length(file_name)
    
    % Displacement from tracked particles on deformed frame
    disp_A2B_parCoordB = -uvw_B2A_prev{ImgSeqNum-1};
    parCoordB = parCoord_prev{ImgSeqNum};

    %%%%% Plot displacements %%%%%
    clf, plotCone3(parCoordB(:,1)*xstep,parCoordB(:,2)*ystep,parCoordB(:,3)*zstep, ...
        disp_A2B_parCoordB(:,1)*xstep ,disp_A2B_parCoordB(:,2)*ystep ,disp_A2B_parCoordB(:,3)*zstep  );
    set(gca,'fontsize',18); view(3); box on; axis equal; axis tight;  
    title(['Tracked cumulative displacement (#',num2str(ImgSeqNum),')'],'fontweight','normal');
    xlabel('x'); ylabel('y'); zlabel('z');
    axis([xstep*MPTPara.gridxyzROIRange.gridx(1), xstep*MPTPara.gridxyzROIRange.gridx(2), ...
          ystep*MPTPara.gridxyzROIRange.gridy(1), ystep*MPTPara.gridxyzROIRange.gridy(2), ...
          zstep*MPTPara.gridxyzROIRange.gridz(1), zstep*MPTPara.gridxyzROIRange.gridz(2)]);
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    
end
close(v);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%% Compute trajectory %%%%%

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
            parCoordATraj{parInd}(ImgSeqNum-1,1:3) = parCoordB(track_A2B(parInd),1:3);
        else
            parCoordATraj{parInd}(ImgSeqNum-1,1:3) = [nan,nan,nan]; 
        end
    end
     
end
 

%%%%% Plot tracked trajectories %%%%% 
disp('%%%%% Plot tracked trajectories %%%%%'); fprintf('\n');
figure,
for parInd = 1:size(parCoordA,1)
    try
        wayPoints = parCoordATraj{parInd}; 
        if (size(resultDisp,1)+1)<4
            hold on; line(wayPoints(isnan(wayPoints(:,1))<1,1),wayPoints(isnan(wayPoints(:,1))<1,2),wayPoints(isnan(wayPoints(:,1))<1,3)); view(3); % straight lines
        else
            hold on; fnplt(cscvn(wayPoints(isnan(wayPoints(:,1))<1,:)'),'',1);
        end
        % if sum(1-isnan(wayPoints(:,1)))>1 % Don't show if there is only one point on the trajectory
        %    hold on; plot3(wayPoints(:,1),wayPoints(:,2),wayPoints(:,3),'.','markersize',8);
        % end
    catch
    end
end

set(gca,'fontsize',18); view(3); box on; axis equal; axis tight;  
title('Tracked particle trajectory','fontweight','normal');
xlabel('x'); ylabel('y'); zlabel('z');
axis([xstep*MPTPara.gridxyzROIRange.gridx(1), xstep*MPTPara.gridxyzROIRange.gridx(2), ...
      ystep*MPTPara.gridxyzROIRange.gridy(1), ystep*MPTPara.gridxyzROIRange.gridy(2), ...
      zstep*MPTPara.gridxyzROIRange.gridz(1), zstep*MPTPara.gridxyzROIRange.gridz(2)]);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Press "Ctrl + C" to modify codes below to plot interpolated ...' ...
      'displacements and strains on a uniform grid mesh']);
disp(['Press "Enter" key to keep running the code']);
pause; 


ImgSeqNum = 2; % TODO: assign a Frame #
 

%%%%% Previously tracked displacement field %%%%%
resultDispCurr = resultDisp{ImgSeqNum-1};
resultDefGradCurr = resultDefGrad{ImgSeqNum-1};
disp_A2B_parCoordB = resultDispCurr.disp_A2B_parCoordB;
parCoordB = resultDispCurr.parCoordB;
  
%%%%% Shift rigid body translations %%%%%
disp_A2B_parCoordB(:,1) = disp_A2B_parCoordB(:,1) - median(disp_A2B_parCoordB(:,1));
disp_A2B_parCoordB(:,2) = disp_A2B_parCoordB(:,2) - median(disp_A2B_parCoordB(:,2));
disp_A2B_parCoordB(:,3) = disp_A2B_parCoordB(:,3) - median(disp_A2B_parCoordB(:,3));


%%%%% Interpolate scatterred data to gridded data %%%%%
addpath('./Scatter2Grid3D/');
sxyz = min([round(0.5*MPTPara.f_o_s),20])*[1,1,1]; % Step size for griddata
smoothness = 1e-3; % Smoothness for regularization; "smoothness=0" means no regularization

[x_Grid_refB,y_Grid_refB,z_Grid_refB,u_Grid_refB]=funScatter2Grid3D(parCoordB(:,1),parCoordB(:,2),parCoordB(:,3),disp_A2B_parCoordB(:,1),sxyz,smoothness);
[~,~,~,v_Grid_refB]=funScatter2Grid3D(parCoordB(:,1),parCoordB(:,2),parCoordB(:,3),disp_A2B_parCoordB(:,2),sxyz,smoothness);
[~,~,~,w_Grid_refB]=funScatter2Grid3D(parCoordB(:,1),parCoordB(:,2),parCoordB(:,3),disp_A2B_parCoordB(:,3),sxyz,smoothness);
 
% Build a displacement vector
uvw_Grid_refB_Vector=[u_Grid_refB(:),v_Grid_refB(:),w_Grid_refB(:)]'; uvw_Grid_refB_Vector=uvw_Grid_refB_Vector(:);

% Calculate deformation gradient
D_Grid = funDerivativeOp3(size(x_Grid_refB,1),size(x_Grid_refB,2),size(x_Grid_refB,3),sxyz); % Central finite difference operator
F_Grid_refB_Vector=D_Grid*uvw_Grid_refB_Vector; % {F}={D}{U}

% Change "uvw_Grid_refB_Vector" and "F_Grid_refB_Vector" in physical world units
uvw_Grid_refB_Vector_PhysWorld = [u_Grid_refB(:)*xstep,v_Grid_refB(:)*ystep,w_Grid_refB(:)*zstep]';
uvw_Grid_refB_Vector_PhysWorld = uvw_Grid_refB_Vector_PhysWorld(:);

F_Grid_refB_Vector_PhysWorld = [F_Grid_refB_Vector(1:9:end), ...             % F11
                                F_Grid_refB_Vector(2:9:end)*ystep/xstep, ... % F21
                                F_Grid_refB_Vector(3:9:end)*zstep/xstep, ... % F31
                                F_Grid_refB_Vector(4:9:end)*xstep/ystep, ... % F12
                                F_Grid_refB_Vector(5:9:end), ...             % F22
                                F_Grid_refB_Vector(6:9:end)*zstep/ystep, ... % F32
                                F_Grid_refB_Vector(7:9:end)*xstep/zstep, ... % F13
                                F_Grid_refB_Vector(8:9:end)*ystep/zstep, ... % F23
                                F_Grid_refB_Vector(9:9:end)]';               % F33
F_Grid_refB_Vector_PhysWorld = F_Grid_refB_Vector_PhysWorld(:);

%%%%% Cone plot grid data: displecement %%%%%
figure, plotCone3(x_Grid_refB*xstep,y_Grid_refB*ystep,z_Grid_refB*zstep,u_Grid_refB*xstep,v_Grid_refB*ystep,w_Grid_refB*zstep);
set(gca,'fontsize',18); view(3); box on; axis equal; axis tight;  
title('Tracked cumulative displacement','fontweight','normal');
axis([xstep*MPTPara.gridxyzROIRange.gridx(1), xstep*MPTPara.gridxyzROIRange.gridx(2), ...
      ystep*MPTPara.gridxyzROIRange.gridy(1), ystep*MPTPara.gridxyzROIRange.gridy(2), ...
      zstep*MPTPara.gridxyzROIRange.gridz(1), zstep*MPTPara.gridxyzROIRange.gridz(2)]);

%%%%% Generate an FE-mesh %%%%%
[coordinatesFEM_refB,elementsFEM_refB] = funMeshSetUp3(x_Grid_refB*xstep,y_Grid_refB*ystep,z_Grid_refB*zstep);

%%%%% Cone plot grid data: displacement %%%%%
Plotdisp_show3(uvw_Grid_refB_Vector_PhysWorld, coordinatesFEM_refB*diag([xstep,ystep,zstep]), elementsFEM_refB,[],'NoEdgeColor');
 
%%%%% Cone plot grid data: infinitesimal strain %%%%%
Plotstrain_show3(F_Grid_refB_Vector_PhysWorld, coordinatesFEM_refB*diag([xstep,ystep,zstep]), elementsFEM_refB,[],'NoEdgeColor',xstep,tstep);
 

 










