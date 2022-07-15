% %%%%%%%%%%%%%%%%%% SerialTrack (3D incremental mode) %%%%%%%%%%%%%%%%%
% Main file of code "SerialTrack"
% ***********************************************
% Dimension:            3D
% Particle rigidity:    hard 
% Tracking mode:        incremental
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

ImgSeqNum=1; [file_name,Img] = funReadImage3(fileNameAll,ImgSeqNum); % Load reference image

try if isempty(fileFolder)~=1, cd(SerialTrackPath); end; catch; end % Come back to the main path

%%%%% Update MPTPara %%%%%
MPTPara.gridxyzROIRange.gridx = [1,size(Img{1},1)];
MPTPara.gridxyzROIRange.gridy = [1,size(Img{1},2)];
MPTPara.gridxyzROIRange.gridz = [1,size(Img{1},3)];

% figure, imagesc3(Img{1}) % To check/plot volumetric image
disp('%%%%%% Load reference image: Done! %%%%%%'); fprintf('\n');

%%%%% Load image mask file %%%%%
try load(im_roi_mask_file_path); catch; end
try MPTPara.ImgRefMask = im_roi'; % Load stored image roi if existed
catch, MPTPara.ImgRefMask = ones(size(Img{1})); % Set up default image mask file
end
disp('%%%%%% Load image mask file: Done! %%%%%%'); fprintf('\n');
  
   
%% ====== Detect and localize particles ======
%%%%% Particle detection and localization parameters %%%%%
%%%%% Bead Parameter %%%%%
% BeadPara.detectionMethod = 2;   % Particle detection method: 1 = TPT (blob finding + radial projection), 
%                                                              2 = TracTrac (LoG blob finding + lsq fit of gaussian)
% BeadPara.thres = 0.4;           % Threshold for detecting particles
% BeadPara.beadRad = 0;           % Estimated radius of a single particle [px]
% BeadPara.minSize = 2;           % Minimum volume of a single particle [px^3]
% BeadPara.maxSize = 1000;        % Maximum volume of a single particle [px^3]
% BeadPara.winSize = [5, 5, 5];   % Default [not currently used]
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
    % x{1}{ImgSeqNum} = radialcenter3dvec(double(Img{ImgSeqNum}),x{1}{ImgSeqNum},BeadPara); % Localize particles
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Store particle positions as "parCoordA" %%%%%
x{1}{ImgSeqNum} = x{1}{ImgSeqNum} + [MPTPara.gridxyzROIRange.gridx(1)-1, MPTPara.gridxyzROIRange.gridy(1)-1, MPTPara.gridxyzROIRange.gridz(1)-1];
parCoordA = x{1}{ImgSeqNum};

%%%%% Remove bad parCoord outside the image area %%%%%
for tempi=1:3, parCoordA( parCoordA(:,tempi)>size(Img{ImgSeqNum},tempi), : ) = []; end
for tempi=1:3, parCoordA( parCoordA(:,tempi)<1, : ) = []; end
 
%%%%% Plot %%%%%
figure, plot3(parCoordA(:,1),parCoordA(:,2),parCoordA(:,3),'bo');
view(3); box on; axis equal; axis tight; set(gca,'fontsize',18); 
title('Detected particles in ref image','fontweight','normal');
 
%%%%% Report detected beads # %%%%%
disp(['Detected particle # in ref image: ',num2str(size(parCoordA,1))]);
disp('%%%%%% Detect particles: Done! %%%%%%'); fprintf('\n');
 


%% %%%%% Initialization %%%%%
%%%%% MPT Parameter %%%%%
% MPTPara.f_o_s = 60;              % Size of search field: max(|u|,|v|,|w|)
% MPTPara.n_neighborsMax = 25;     % Max # of neighboring particles
% MPTPara.n_neighborsMin = 1;      % Min # of neighboring particles
% MPTPara.gbSolver = 2;            % Global step solver: 1-moving least square fitting; 2-global regularization; 3-ADMM iterations
% MPTPara.smoothness = 1e-1;       % Coefficient of regularization
% MPTPara.outlrThres = 5;          % Threshold for removing outliers in TPT
% MPTPara.maxIterNum = 20;         % Max ADMM iteration number
% MPTPara.iterStopThres = 1e-2;    % ADMM iteration stopping threshold
% MPTPara.strain_n_neighbors = 20; % # of neighboring particles used in strain gauge
% MPTPara.strain_f_o_s = 60;       % Size of virtual strain gauge
% MPTPara.usePrevResults = 0;      % Whether use previous results or not: 0-no; 1-yes;
% MPTPara.distMissing = 5;         % Distance threshold to check whether particle has a match or not [px]

%%%%%% To store results %%%%%
parCoord_prev = cell(length(file_name)-1,1);    parCoord_prev{1} = parCoordA;
track_A2B_prev = cell(length(file_name)-1,1);   track_B2A_prev = cell(length(file_name)-1,1);
uvw_B2A_prev = cell(length(file_name)-1,1);

 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ImgSeqNum = 2 : length(file_name)  % "ImgSeqNum" is the frame index
    
    disp(['====== Frame #',num2str(ImgSeqNum),' ======']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Load image volumetric data %%%%%
    try if isempty(fileFolder)~=1, cd(fileFolder); end; catch; end % Open image folder
    tempvol = load(file_name{ImgSeqNum}); fieldName = fieldnames(tempvol);
    Img{2} = getfield(tempvol,fieldName{1}); clear tempvol; %#ok<GFLD>
    if iscell(Img{2}), Img{2}=Img{2}{1}; end
    try if isempty(fileFolder)~=1, cd(SerialTrackPath); end; catch; end % Come back to the main path

    %%%%% SerialTrack particle tracking %%%%%
    [parCoordB_temp,uvw_B2A_temp,~,~,track_A2B_temp,track_B2A_temp] = fun_SerialTrack_3D_HardPar( ...
       ImgSeqNum,Img{2},BeadPara,MPTPara,parCoord_prev{ImgSeqNum-1},parCoord_prev(2:end),uvw_B2A_prev);
     
    %%%%% Store results %%%%%
    parCoord_prev{ImgSeqNum} = parCoordB_temp;
    uvw_B2A_prev{ImgSeqNum-1} = uvw_B2A_temp;  % incremental displacement
    track_A2B_prev{ImgSeqNum-1} = track_A2B_temp;
    track_B2A_prev{ImgSeqNum-1} = track_B2A_temp;
    
end
  

%% %%%%% Incremental tracking ratio %%%%%
disp('%%%%% Calculate incremental tracking ratio %%%%%'); fprintf('\n');
track_ratio = zeros(length(file_name)-1,1);
DefType = 'exp'; defList = [2:1:length(file_name)]';
  
for ImgSeqNum = 2 : length(file_name)
    track_A2B = track_A2B_prev{ImgSeqNum-1}; 
    track_ratio(ImgSeqNum-1) = length(track_A2B(track_A2B>0))/size(parCoord_prev{ImgSeqNum},1);      
end
 
fig=figure; ax=axes; hold on; plot(defList,track_ratio,'r^-.','linewidth',1);
adjust_fig(fig,ax,'','',''); box on; title('');
xlabel('Frame #'); ylabel('Tracking ratio');
try axis([2,length(file_name),0,1]); catch, end

%%%%% Save results %%%%%
disp('%%%%%% SerialTrack 3D hard particle tracking: Done! %%%%%%'); fprintf('\n');
results_file_name = 'results_3D_hardpar.mat';
mkdir results
save(['./results/' results_file_name],'parCoord_prev','uvw_B2A_prev','track_A2B_prev','track_B2A_prev');
 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%% Visualize tracked incremental displacement of each frame %%%%%
disp('%%%%% Plot tracked incremental deformations %%%%%'); fprintf('\n');

%%%%% Experimental parameters %%%%%
try xstep = MPTPara.xstep; catch, xstep = 1; end % unit: um/px
try ystep = MPTPara.ystep; catch, ystep = xstep; end
try zstep = MPTPara.zstep; catch, zstep = xstep; end
try tstep = MPTPara.tstep; catch, tstep = 1; end % unit: us  
 
%%%%% Plot tracked incremental displacement field %%%%%
%%%%% Make a video %%%%%
v = VideoWriter('video_3D_inc.avi'); v.FrameRate = 5; open(v); figure,

for ImgSeqNum = 2:length(file_name) % ImgSeqNum: Frame #
    
    % Displacement from tracked particles on deformed frame
    disp_A2B_parCoordB = -uvw_B2A_prev{ImgSeqNum-1};
    parCoordB = parCoord_prev{ImgSeqNum};

    %%%%% Plot displacements %%%%%
    clf, plotCone3(parCoordB(:,1)*xstep, parCoordB(:,2)*ystep, parCoordB(:,3)*zstep, ...
        disp_A2B_parCoordB(:,1)*xstep/tstep, disp_A2B_parCoordB(:,2)*ystep/tstep, disp_A2B_parCoordB(:,3)*zstep/tstep );
    set(gca,'fontsize',18); box on; axis equal; view(3);
    title(['Tracked velocity (#',num2str(ImgSeqNum),')'],'fontweight','normal');
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
if size(parCoord_prev,1) > 1  % Only for multiple frames

disp('%%%%% Merge tracked trajectory segments %%%%%'); fprintf('\n');
 
%%%%% Initialization %%%%%
try xstep = MPTPara.xstep; catch, xstep = 1; end      % um2px ratio
try ystep = MPTPara.ystep; catch, ystep = xstep; end
try zstep = MPTPara.zstep; catch, zstep = xstep; end
try tstep = MPTPara.tstep; catch tstep = 1; end % time gap between consecutive frames
trajInd = 0; % index of trajectory segments

%%%%% These parameters are assigned values in "Example_main_xxx.m" file
% distThres = 3;           % distance threshold to connect split trajectory segments
% extrapMethod = 'pchip';  % extrapolation scheme to connect split trajectory segments
%                          % suggestion: 'nearest' for Brownian motion                          
% minTrajSegLength = 4;    % the minimum length of trajectory segment that will be extrapolate 
% maxGapTrajSeqLength = 0; % the max frame# gap between connected trajectory segments


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Compute and collect all trajectory segments %%%%%
for tempk = 1 : size(parCoord_prev,1)  % Find trajectories passing through particles in frame #tempk
    
    parCoordCurr = parCoord_prev{tempk}; % Find trajectories passing parCoordCurr
    clear parCoordTrajCurr; parCoordTrajCurr = cell(size(parCoordCurr,1),1); % Initialize a cell structure to store #tempk trajectories
    
    % Add particles in frame #tempk to "parCoordTrajCurr"
    for tempj = 1:length(parCoordCurr)
        parCoordTrajCurr{tempj}(tempk,1:3) = parCoordCurr(tempj,1:3);
    end
    
    for ImgSeqNum = (1+tempk) : size(parCoord_prev,1) % For each later tracked incremental deformation
        
        parCoordB = parCoord_prev{ImgSeqNum}; % Current particle coordinates in frame #ImgSeqNum
        parCoordB_Ind = []; parCoordCurr_Ind = []; % Initialize particle index
        track_B2A_Curr = track_B2A_prev{ImgSeqNum-1}; % Current tracking in frame #ImgSeqNum
        for tempi = 1:length(track_B2A_Curr)
            try
                % Back propogation
                if ImgSeqNum > 2
                    for tempj = (ImgSeqNum-2) : -1 : tempk
                        track_B2A_Last = track_B2A_prev{tempj};
                        track_B2A_Curr(tempi) = track_B2A_Last(track_B2A_Curr(tempi));
                    end
                end
                % Check whether the trajectory has already been plotted previously
                if tempk>1
                    track_B2A_SecLast = track_B2A_prev{tempk-1};
                    track_B2A_SecLast_temp = track_B2A_SecLast(track_B2A_Curr(tempi)); % if this is not 0, means it's already been plotted
                else % Trajectories from first two frames (tempk=1) will be added
                    track_B2A_SecLast_temp = 0;
                end
                % Assign index values
                if (track_B2A_Curr(tempi)>0) && (track_B2A_SecLast_temp==0)
                    parCoordB_Ind = [parCoordB_Ind; tempi];
                    parCoordCurr_Ind = [parCoordCurr_Ind; track_B2A_Curr(tempi)];
                end
            catch
            end
        end
        
        for tempj = 1:length(parCoordCurr_Ind) % Add found trajectories to cell structure "parCoordTraj"
            parCoordTrajCurr{parCoordCurr_Ind(tempj)}(ImgSeqNum,1:3) = parCoordB(parCoordB_Ind(tempj),1:3);
        end
    end
    
    for parInd = 1:size(parCoordTrajCurr,1)
        wayPoints = parCoordTrajCurr{parInd};
        if ~isempty(wayPoints)
            wayPoints(wayPoints(:,1)==0,:) = wayPoints(wayPoints(:,1)==0,:)*nan;
            wayPoints = [wayPoints; nan(size(parCoord_prev,1)-size(wayPoints,1),3)];
            trajInd = trajInd + 1;
            parCoordTraj{trajInd} = wayPoints;
        end
    end
      
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%% Merge trajectory segments %%%%%
% Find the starting point and length of each trajectory segment
% parCoordTrajPara has size [xxx, 2], and each row is 
% parCoordTrajPara(xxx, 1:2) = [first non-NAN position, trajectory segment length]
parCoordTrajPara = []; parCoordTraj = parCoordTraj(:);
for tempi = 1 : size(parCoordTraj,1)
    parCoordTrajCurr = parCoordTraj{tempi};
    parCoordTrajPara(tempi,1) = find(isnan(parCoordTrajCurr(:,1))==0, 1, 'first');
    parCoordTrajPara(tempi,2) = sum(1-isnan(parCoordTrajCurr(:,1)));
end

%%%%% Try to merge trajectory segments by extrapolating the particle position %%%%%
hbar = waitbar(0,'wait');
for tempMergeTime = 1:4 % Try to merge four times
    
    for tempm = 0:maxGapTrajSeqLength % tempm is the # of missing particles between trajectory segments
        
        for tempk = (size(parCoord_prev,1)-1) : -1 : minTrajSegLength % For trajectory segments with length [(size(parCoord_prev,1)-1) : -1 : minTrajSegLength]
            
            [row,col] = find( parCoordTrajPara(:,2)==tempk ); % Find trajectory segments whose "length==tempk"
            [row1,~] = find( parCoordTrajPara(:,2)<size(parCoord_prev,1)+1-tempk & parCoordTrajPara(:,2)>0 ); % Find trajectory setments with length requirement
             
            for tempi = 1:length(row) % For each trajectory segment whose "length==tempk"
                 
                tempWatibar = tempi/length(row)/4/(maxGapTrajSeqLength+1)/((size(parCoord_prev,1)-1)-minTrajSegLength+1) + ...
                    ((size(parCoord_prev,1)-1)-tempk+1)/4/(maxGapTrajSeqLength+1)/((size(parCoord_prev,1)-1)-minTrajSegLength+1) + ...
                    (tempm)/4/(maxGapTrajSeqLength+1) + (tempMergeTime-1)/4;
                
                waitbar(tempWatibar);
                 
                parCoordTrajMat = cell2mat( parCoordTraj );
                parCoordTrajCurr = parCoordTraj{row(tempi)}; % For each trajectory segment whose "length==tempk"
                parCoordTrajCurrPred_x = fillmissing(parCoordTrajCurr(:,1),extrapMethod); % fill missing using 'pchip' method
                parCoordTrajCurrPred_y = fillmissing(parCoordTrajCurr(:,2),extrapMethod); % fill missing using 'pchip' method
                parCoordTrajCurrPred_z = fillmissing(parCoordTrajCurr(:,3),extrapMethod); % fill missing using 'pchip' method
                % figure, plot3(parCoordTrajCurrPred_x,parCoordTrajCurrPred_y,parCoordTrajCurrPred_z,'-o');
                
                
                %%%%% Find all probable trajectory segments in the positive direction %%%%%
                if (sum(parCoordTrajPara(row(tempi),1:2))+tempm<(size(parCoord_prev,1)+1)) && (sum(parCoordTrajPara(row(tempi),1:2))+tempm>0)
                    
                    [row2,col2] = find( parCoordTrajPara(:,1) == sum(parCoordTrajPara(row(tempi),1:2))+tempm ); % starting point requirement
                    row3 = intersect(row1,row2);
                    temp1 = size(parCoord_prev,1)*(row3-1) +(sum(parCoordTrajPara(row(tempi),1:2)))+tempm; temp1=temp1'; temp1=temp1(:); % Index in parCoordTrajMat
                    temp2 = parCoordTrajMat(temp1, 1:3);
                    
                    % hold on; plot3(temp2(:,1),temp2(:,2),temp2(:,2),'.'); pause;
                    
                    temp3 = sqrt((temp2(:,1)-parCoordTrajCurrPred_x( sum(parCoordTrajPara(row(tempi),1:2))+tempm )).^2 + ...
                                (temp2(:,2)-parCoordTrajCurrPred_y( sum(parCoordTrajPara(row(tempi),1:2))+tempm )).^2 + ...
                                (temp2(:,3)-parCoordTrajCurrPred_z( sum(parCoordTrajPara(row(tempi),1:2))+tempm )).^2 );
                    
                    [temp3min,temp3minind] = min(temp3);
                    if temp3min < distThres % Find the continuing trajectory segment %JY!!!! threshold distance 3
                        
                        % Merge trajectory segment
                        parCoordTraj{row(tempi)}( parCoordTrajPara(row3(temp3minind),1) : parCoordTrajPara(row3(temp3minind),1)+parCoordTrajPara(row3(temp3minind),2)-1, 1:3 ) = ...
                            parCoordTraj{row3(temp3minind)}( parCoordTrajPara(row3(temp3minind),1) : parCoordTrajPara(row3(temp3minind),1)+parCoordTrajPara(row3(temp3minind),2)-1, 1:3 );
                        % Remove repeated trajectory segment
                        parCoordTraj{row3(temp3minind)} = parCoordTraj{row3(temp3minind)}*nan;
                        
                        % Update varaible "parCoordTrajPara" for trajectory segment length
                        parCoordTrajPara(row(tempi),2) = parCoordTrajPara(row(tempi),2) + parCoordTrajPara(row3(temp3minind),2) + tempm;
                        parCoordTrajPara(row3(temp3minind),1:2) = [0,0];
                        
                        %Fillmissing parCoordTraj{row(tempi)}
                        temp = parCoordTraj{row(tempi)}(parCoordTrajPara(row(tempi),1) : sum(parCoordTrajPara(row(tempi),1:2))-1, 1:3);
                        temp_x = fillmissing(temp(:,1),extrapMethod); % fill missing
                        temp_y = fillmissing(temp(:,2),extrapMethod); % fill missing
                        temp_z = fillmissing(temp(:,3),extrapMethod); % fill missing
                        
                        parCoordTraj{row(tempi)}(parCoordTrajPara(row(tempi),1) : sum(parCoordTrajPara(row(tempi),1:2))-1, 1:3) = [temp_x, temp_y, temp_z];
                        
                    end
                    
                end
                
                
                %%%%% Find all probable trajectory segments in the negative direction %%%%%
                if (sum(parCoordTrajPara(row(tempi),1:2))+tempm<(size(parCoord_prev,1)+1)) && ((parCoordTrajPara(row(tempi),1)-1)-tempm>0)
                    
                    [row2,col2] = find( sum(parCoordTrajPara,2) == parCoordTrajPara(row(tempi),1)-tempm ); % ending point requirement
                    row3 = intersect(row1,row2);
                    temp1 = size(parCoord_prev,1)*(row3-1) + (parCoordTrajPara(row(tempi),1)-1)-tempm; temp1=temp1'; temp1=temp1(:); % Index in parCoordTrajMat
                    temp2 = parCoordTrajMat(temp1, 1:3);
                    
                    % hold on; plot(tempyy(:,1),tempyy(:,2),'.');
                    temp3 = sqrt( ( temp2(:,1)-parCoordTrajCurrPred_x( (parCoordTrajPara(row(tempi),1)-1)-tempm ) ).^2 + ...
                                ( temp2(:,2)-parCoordTrajCurrPred_y( (parCoordTrajPara(row(tempi),1)-1)-tempm ) ).^2 + ...
                                ( temp2(:,3)-parCoordTrajCurrPred_z( (parCoordTrajPara(row(tempi),1)-1)-tempm ) ).^2  );
                    
                    [temp3min,temp3minind] = min(temp3);
                    if temp3min < distThres % Find the continuing trajectory segment  
                        % Merge trajectory segment
                        parCoordTraj{row(tempi)}( parCoordTrajPara(row3(temp3minind),1) : parCoordTrajPara(row3(temp3minind),1)+parCoordTrajPara(row3(temp3minind),2)-1, 1:3 ) = ...
                            parCoordTraj{row3(temp3minind)}( parCoordTrajPara(row3(temp3minind),1) : parCoordTrajPara(row3(temp3minind),1)+parCoordTrajPara(row3(temp3minind),2)-1, 1:3 );
                        % Remove repeated trajectory segment
                        parCoordTraj{row3(temp3minind)} = parCoordTraj{row3(temp3minind)}*nan;
                        % Update varaible "parCoordTrajPara" for both trajectory segment starting point and its length
                        parCoordTrajPara(row(tempi),2) = parCoordTrajPara(row(tempi),2) + parCoordTrajPara(row3(temp3minind),2) + tempm;
                        parCoordTrajPara(row(tempi),1) = parCoordTrajPara(row3(temp3minind),1);
                        parCoordTrajPara(row3(temp3minind),1:2) = [0,0];
                        
                        %Fillmissing parCoordTraj{row(tempi)}
                        temp = parCoordTraj{row(tempi)}(parCoordTrajPara(row(tempi),1) : sum(parCoordTrajPara(row(tempi),1:2))-1, 1:3);
                        temp_x = fillmissing(temp(:,1),extrapMethod); % fill missing
                        temp_y = fillmissing(temp(:,2),extrapMethod); % fill missing
                        temp_z = fillmissing(temp(:,3),extrapMethod); % fill missing
                        parCoordTraj{row(tempi)}(parCoordTrajPara(row(tempi),1) : sum(parCoordTrajPara(row(tempi),1:2))-1, 1:3) = [temp_x, temp_y, temp_z];
                        
                    end 
                    
                end
                
            end % End of "for each trajectory segment whose "length==tempk" "
            
            
            
        end % End of "for trajectory segments with length [(size(parCoord_prev,1)-1) : -1 : minTrajSegLength]"
        
    end % End of tempm
    
end % End of "for tempMergeTime = 1:5 % Try to merge four times"
close(hbar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plot tracked trajectories %%%%%
disp('%%%%% Plot tracked trajectories %%%%%'); fprintf('\n');

figure,
% CN = [36 63 146; 20 76 156; 28 100 175; 9 115 186; 23 128 196; 33 148 210;
%     25 158 218; 19 172 226; 25 186 229; 69 196 221; 118 205 214; 157 216 215;
%     169 220 217; 193 229 224; 216 237 233; 234 246 245]/255;
% CN = CN(end-size(parCoord_prev,1):end, 1:3); % Colormap I used in the manuscript


for tempi = 1:size(parCoordTrajPara,1)
    
    wayPoints = parCoordTraj{tempi};
    
    if (length(wayPoints(isnan(wayPoints(:,1))<1,1))+1)<4
        hold on; line(xstep*wayPoints(isnan(wayPoints(:,1))<1,1), ...
                      ystep*wayPoints(isnan(wayPoints(:,1))<1,2), ...
                      zstep*wayPoints(isnan(wayPoints(:,1))<1,3), 'linewidth', 1); view(2); % straight lines
    else
        hold on; fnplt(cscvn(diag([xstep,ystep,zstep])*wayPoints(isnan(wayPoints(:,1))<1,:)'),'',1);
    end
     
    %%%%% Codes to plot trajectories with frame-dependent color %%%%% 
    % if sum(1-isnan(wayPoints(:,1)))>1  % Don't show if there is only one point on the trajectory
    %     hold on; plot3(xstep*wayPoints(:,1),ystep*wayPoints(:,2),zstep*wayPoints(:,3),'r.','markersize',5);
    % end
    % 
    % for tempj = 1:size(parCoord_prev,1)-1
    %     hold on; line(xstep*[wayPoints(tempj,1),wayPoints(tempj+1,1)], ...
    %                   ystep*[wayPoints(tempj,2),wayPoints(tempj+1,2)], ...
    %                   zstep*[wayPoints(tempj,3),wayPoints(tempj+1,3)], 'linewidth',1.2, 'color', CN(tempj,:) );
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
end
 
set(gca,'fontsize',18); view(3); box on; axis equal; axis tight;  
title('Tracked particle trajectory','fontweight','normal');
xlabel('x'); ylabel('y'); zlabel('z');
axis([xstep*MPTPara.gridxyzROIRange.gridx(1), xstep*MPTPara.gridxyzROIRange.gridx(2), ...
      ystep*MPTPara.gridxyzROIRange.gridy(1), ystep*MPTPara.gridxyzROIRange.gridy(2), ...
      zstep*MPTPara.gridxyzROIRange.gridz(1), zstep*MPTPara.gridxyzROIRange.gridz(2)]);

  
%% %%%%% Compute cumulative tracking ratio from emerged trajectories %%%%%
 
disp('%%%%% Plot tracked umulative displacements %%%%%'); fprintf('\n');
parCoordTrajMat = cell2mat( parCoordTraj );

[row1,col1] = find(isnan(parCoordTrajMat(1:length(file_name):end,1))==0);
trackParaccum_ind = row1;
trackParaccum_track_ratio = [];

for ImgSeqNum = 2:length(file_name)
    [row2,col2] = find(isnan(parCoordTrajMat(ImgSeqNum:length(file_name):end,1))==0);
    trackParaccum_ind = intersect(row2,trackParaccum_ind);
    trackParaccum_track_ratio(ImgSeqNum-1) = length(trackParaccum_ind) / size(parCoord_prev{1},1);
end

defList = [2:1:length(file_name)];
fig=figure; ax=axes; hold on; plot(defList,trackParaccum_track_ratio,'bs--','linewidth',1);
adjust_fig(fig,ax,'','',''); box on; title('');
xlabel('Image #'); ylabel('Tracking ratio');
try axis([2,length(file_name),0,1]); catch; end

%%%%% Plot tracked cumulative displacement field %%%%%
%%%%% Make a video %%%%%
v = VideoWriter('video_3D_inc_accum.avi'); v.FrameRate = 5; open(v); figure,
for ImgSeqNum = 2:length(file_name)
    
    parCoordA = parCoordTrajMat(1:length(file_name):end,1:3);
    parCoordB = parCoordTrajMat(ImgSeqNum:length(file_name):end,1:3);
    parCoordAaccum = parCoordA(trackParaccum_ind,:);
    parCoordBaccum = parCoordB(trackParaccum_ind,:);
    disp_A2Baccum = parCoordBaccum - parCoordAaccum;
    
    % ----- Cone plot grid data: displecement -----
    clf; plotCone3(xstep*parCoordBaccum(:,1),ystep*parCoordBaccum(:,2),zstep*parCoordBaccum(:,3), ...
                   xstep*disp_A2Baccum(:,1),ystep*disp_A2Baccum(:,2),zstep*disp_A2Baccum(:,3));
    set(gca,'fontsize',18); view(3); box on; axis equal; axis tight; 
    title(['Tracked cumulative disp (#',num2str(ImgSeqNum),')'],'fontweight','normal');
    xlabel('x'); ylabel('y'); zlabel('z');
    axis([xstep*MPTPara.gridxyzROIRange.gridx(1), xstep*MPTPara.gridxyzROIRange.gridx(2), ...
          ystep*MPTPara.gridxyzROIRange.gridy(1), ystep*MPTPara.gridxyzROIRange.gridy(2), ...
          zstep*MPTPara.gridxyzROIRange.gridz(1), zstep*MPTPara.gridxyzROIRange.gridz(2)]);
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    
end
close(v);

end    % End of if:  Only work for multiple frames


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Press "Ctrl + C" to modify codes below to plot interpolated ...' ...
      'displacements and strains on a uniform grid mesh']);
disp(['Press "Enter" key to keep running the code']);
pause; 


ImgSeqNum = 2; % TODO: assign a Frame #


%%%%% Previously tracked displacement field %%%%%
if size(parCoord_prev,1) == 1 % if there are just two frames
    parCoordBaccum = parCoordB;
    parCoordAaccum = parCoordB - disp_A2B_parCoordB;
    disp_A2Baccum = parCoordBaccum - parCoordAaccum; 
else % for multiple frames
    parCoordA = parCoordTrajMat(1:length(file_name):end,1:3);
    parCoordB = parCoordTrajMat(ImgSeqNum:length(file_name):end,1:3);
    parCoordAaccum = parCoordA(trackParaccum_ind,1:3);
    parCoordBaccum = parCoordB(trackParaccum_ind,1:3);
    disp_A2Baccum = parCoordBaccum - parCoordAaccum;
end

%%%%% Interpolate scatterred data to gridded data %%%%%
addpath('./Scatter2Grid3D/'); % MATLAB Exchange File (see Ref[5])
sxyz = min([round(0.5*MPTPara.f_o_s),20])*[1,1,1]; % Step size for griddata
smoothness = 1e-3; % Smoothness for regularization; "smoothness=0" means no regularization

[x_Grid_refB,y_Grid_refB,z_Grid_refB,u_Grid_refB]=funScatter2Grid3D(parCoordBaccum(:,1),parCoordBaccum(:,2),parCoordBaccum(:,3),disp_A2Baccum(:,1),sxyz,smoothness);
[~,~,~,v_Grid_refB]=funScatter2Grid3D(parCoordBaccum(:,1),parCoordBaccum(:,2),parCoordBaccum(:,3),disp_A2Baccum(:,2),sxyz,smoothness);
[~,~,~,w_Grid_refB]=funScatter2Grid3D(parCoordBaccum(:,1),parCoordBaccum(:,2),parCoordBaccum(:,3),disp_A2Baccum(:,3),sxyz,smoothness);

% Apply ROI image mask %% TODO
% [u_Grid_refB, v_Grid_refB] = funRmROIOutside(x_Grid_refB,y_Grid_refB,MPTPara.ImgRefMask,u_Grid_refB,v_Grid_refB);

% Build a displacement vector
uvw_Grid_refB_Vector = [u_Grid_refB(:),v_Grid_refB(:),w_Grid_refB(:)]'; uvw_Grid_refB_Vector=uvw_Grid_refB_Vector(:);

% Calculate deformation gradient
D_Grid = funDerivativeOp3(size(x_Grid_refB,1),size(x_Grid_refB,2),size(x_Grid_refB,3), sxyz); % Central finite difference operator
F_Grid_refB_Vector = D_Grid*uvw_Grid_refB_Vector; % {F}={D}{U}

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
figure, plotCone3(x_Grid_refB*xstep,y_Grid_refB*ystep,z_Grid_refB*zstep,u_Grid_refB*xstep,v_Grid_refB*ystep,w_Grid_refB*zstep );
set(gca,'fontsize',18); view(3); box on; axis equal; axis tight;  
title('Tracked cumulative displacement','fontweight','normal');
axis([xstep*MPTPara.gridxyzROIRange.gridx(1), xstep*MPTPara.gridxyzROIRange.gridx(2), ...
      ystep*MPTPara.gridxyzROIRange.gridy(1), ystep*MPTPara.gridxyzROIRange.gridy(2), ...
      zstep*MPTPara.gridxyzROIRange.gridz(1), zstep*MPTPara.gridxyzROIRange.gridz(2)]);

  
%%%%% Generate an FE-mesh %%%%%
[coordinatesFEM_refB,elementsFEM_refB] = funMeshSetUp3(x_Grid_refB*xstep,y_Grid_refB*ystep,z_Grid_refB*zstep);

%%%%% Cone plot grid data: displacement %%%%%
Plotdisp_show3(uvw_Grid_refB_Vector_PhysWorld, ...
    coordinatesFEM_refB*diag([xstep,ystep,zstep]), elementsFEM_refB,[],'NoEdgeColor');
 
%%%%% Cone plot grid data: infinitesimal strain %%%%%
Plotstrain_show3(F_Grid_refB_Vector_PhysWorld, coordinatesFEM_refB*diag([xstep,ystep,zstep]), ...
    elementsFEM_refB,[],'NoEdgeColor',xstep,tstep);
 

 










