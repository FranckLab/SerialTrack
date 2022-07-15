%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SerialTrack execute main file
% ===================================================
% Dimension:            2D
% Particle rigidity:    hard
% Tracking mode:        incremental
% Deformation mode:     provided particle coordinates 
%
% ===================================================
% Author: Jin Yang, Ph.D.
% Email: jyang526@wisc.edu -or-  aldicdvc@gmail.com
% Date: 07/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
close all; clear all; clc; clearvars -global
disp('************************************************');
disp('*** Welcome to SerialTrack Particle Tracking ***');
disp('************************************************');
addpath( './function/','./src/' );
 

%% Load tracked particle coordinates "parCoordA" and "parCoordB"

fprintf('\n'); disp('************************************************');
disp('This example is to demo how SerialTrack links');
disp('particles whose centroids are already provided.');
disp('************************************************'); fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Option I: If you only have two images %%%%%
% % TODO: load particle centroid locations in the ref image as "parCoordA"
% load('.\Example_main_files\data_lung_flow_parCoordA.mat');
% % TODO: load particle centroid locations in the def image as "parCoordB"
% load('.\Example_main_files\data_lung_flow_parCoordB.mat');
% 
% parCoord_prev{1} = parCoordA;
% parCoord_prev{2} = parCoordB;
% 
% %%%%% Multiple particle tracking (MPT) Parameter %%%%%
% MPTPara.f_o_s = 30;              % Size of search field: max(|u|,|v|,|w|) [px]
% MPTPara.n_neighborsMax = 25;     % Max # of neighboring particles
% MPTPara.n_neighborsMin = 1;      % Min # of neighboring particles
% MPTPara.locSolver = 1;           % Local solver: 1-topology-based feature; 
%                                  %               2-histogram-based feature first and then topology-based feature;
% MPTPara.gbSolver = 3;            % Global step solver: 1-moving least square fitting; 
%                                  %                     2-global regularization; 
%                                  %                     3-ADMM iterations
% MPTPara.smoothness = 1e-2;       % Coefficient of regularization
% MPTPara.outlrThres = 2;          % Threshold for removing outliers in MPT
% MPTPara.maxIterNum = 20;         % Max ADMM iteration number
% MPTPara.iterStopThres = 1e-2;    % ADMM iteration stopping threshold
% MPTPara.strain_n_neighbors = 20; % # of neighboring particles used in strain gauge
% MPTPara.strain_f_o_s = 60;       % Size of virtual strain gauge [px]
% MPTPara.usePrevResults = 0;      % Whether use previous results or not: 0-no; 1-yes;
% MPTPara.distMissing = 2;         % Distance threshold to check whether particle has a match or not



%%%%% Option II: If you have more than two images %%%%%
% TODO: load particle centroid locations and name it as a cell array "parCoord_prev"
load('.\Example_main_files\data_pipe_parCoord.mat');
parCoordA = parCoord_prev{1}; 

% =============================== Please also update MPTPara as:
%%%%%% Multiple particle tracking (MPT) Parameter %%%%%
MPTPara.f_o_s = 15;              % Size of search field: max(|u|,|v|) [px]
MPTPara.n_neighborsMax = 25;     % Max # of neighboring particles
MPTPara.n_neighborsMin = 1;      % Min # of neighboring particles
MPTPara.locSolver = 1;           % Local solver: 1-topology-based feature; 
                                 %               2-histogram-based feature first and then topology-based feature;
MPTPara.gbSolver = 3;            % Global step solver: 1-moving least square fitting; 
                                 %                     2-global regularization; 
                                 %                     3-ADMM iterations
MPTPara.smoothness = 1e-2;       % Coefficient of regularization
MPTPara.outlrThres = 2;          % Threshold for removing outliers in TPT
MPTPara.maxIterNum = 20;         % Max ADMM iteration number
MPTPara.iterStopThres = 1e-2;    % ADMM iteration stopping threshold
MPTPara.strain_n_neighbors = 20; % # of neighboring particles used in strain gauge
MPTPara.strain_f_o_s = 60;       % Size of virtual strain gauge [px]
MPTPara.usePrevResults = 0;      % Whether use previous results or not: 0-no; 1-yes;  
MPTPara.distMissing = 2;         % Distance threshold to check whether particle has a match or not [px]
 

disp('%%%%%% Load tracked particle corodinates: Done! %%%%%%'); fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Image binary mask file %%%%%
MaskFileLoadingMode = 0; % {0: No mask file
                         %  1: Load only one mask file for all frames; 
                         %  2: Load one mask file for each single frame;
                         %  3: Load a MATLAB mat file for all frames;

if MaskFileLoadingMode == 3
    im_roi_mask_file_path = '.\img_par2track_lung\im_roi.mat';  % TODO: Path of the mat file to be used as the mask file
else
    im_roi_mask_file_path = '';  % If there is no mask mat file, leave it as empty;
end
disp('%%%%%% Load image mask file: Done! %%%%%%'); fprintf('\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 1: User defined parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Update MPTPara %%%%%
MPTPara.gridxyROIRange.gridx = [min(floor(parCoordA(:,1))), max(ceil(parCoordA(:,1)))];
MPTPara.gridxyROIRange.gridy = [min(floor(parCoordA(:,2))), max(ceil(parCoordA(:,2)))]; 

%%%%% Problem dimension and units %%%%%
MPTPara.DIM = 2;      % problem dimension
MPTPara.xstep = 1;    % unit: um/px
MPTPara.tstep = 1;    % unit: us -or- 1/frame

%%%%% Code mode %%%%%
MPTPara.mode = 'inc'; % {'inc': incremental mode;
                      %  'accum': cumulative mode}

%%%%% Particle rigidity %%%%%
MPTPara.parType = 'hard'; % {'hard': hard particle;
                          %  'soft': soft particle}

disp('************************************************');
disp(['Dimention: ',num2str(MPTPara.DIM)]);
disp(['Tracking mode: ',MPTPara.mode]);
disp(['Particle type: ',MPTPara.parType]);
disp('************************************************'); fprintf('\n');

%%%% Postprocessing: merge trajectory segments %%%%%
distThres = 1;              % distance threshold to connect split trajectory segments [px]
extrapMethod = 'pchip';     % extrapolation scheme to connect split trajectory segments
% suggestion: 'nearest' for Brownian motion
minTrajSegLength = 10;      % the minimum length of trajectory segment that will be extrapolate [px]
maxGapTrajSeqLength = 0;    % the max frame# gap between connected trajectory segments


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 2: Particle Linking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== Track A & B neighbor-topology match ======
%
%  Directly run this section if coordinates of detected particles are known:
%  Coordinates of particles in reference image:   parCoordA
%  Coordinates of particles in deformed image:    parCoordB
%
%    \  |  /                  \  |  /
%     \ | /                    \ | /
%   --- A ---       v.s.     --- B ---
%      / \                      / \
%     /   \                    /   \
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% To store results %%%%%
track_A2B_prev = cell(length(parCoord_prev)-1,1);
track_B2A_prev = cell(length(parCoord_prev)-1,1);
uv_B2A_prev = cell(length(parCoord_prev)-1,1);

for ImgSeqNum = 2 : length(parCoord_prev) % "ImgSeqNum" is the frame index

    disp(['====== Frame #',num2str(ImgSeqNum),' ======']);

    [matches_A2B,u_B2A_curr_refB,track_A2B] = f_track_serial_match2D( ...
        parCoord_prev{ImgSeqNum-1}, parCoord_prev{ImgSeqNum}, ...
        'f_o_s',MPTPara.f_o_s, 'n_neighborsMax', MPTPara.n_neighborsMax, 'n_neighborsMin', MPTPara.n_neighborsMin, ...
        'locSolver',MPTPara.locSolver,'gbSolver', MPTPara.gbSolver, 'smoothness', MPTPara.smoothness, ...
        'outlrThres', MPTPara.outlrThres, 'maxIterNum', MPTPara.maxIterNum, ...
        'iterStopThres', MPTPara.iterStopThres, 'usePrevResults', MPTPara.usePrevResults, ...
        'strain_n_neighbors',MPTPara.strain_n_neighbors, 'strain_f_o_s',MPTPara.strain_f_o_s, ...
        'gridxyROIRange',MPTPara.gridxyROIRange, 'parCoordB_prev',parCoord_prev(2:end), ...
        'uv_B2A_prev',uv_B2A_prev, 'ImgSeqNum',ImgSeqNum, ...
        'MPTParaDistMissing',MPTPara.distMissing);


    %%%%% Compute track_B2A %%%%%
    matches_A2B = matches_A2B(matches_A2B(:,2)>0,:); % Build matches_A2B_temp
    track_B2A = zeros(size(parCoord_prev{ImgSeqNum},1), 1);
    for tempi = 1:size(matches_A2B)
        track_B2A(matches_A2B(tempi,2)) = matches_A2B(tempi,1);
    end


    %% %%%%% Plotting %%%%%%
    % Compute displacement from tracked particles on deformed frame
    disp_A2B_parCoordB = -u_B2A_curr_refB;
    figure, plotCone2(parCoord_prev{ImgSeqNum}(:,1),parCoord_prev{ImgSeqNum}(:,2), ...
                      disp_A2B_parCoordB(:,1),disp_A2B_parCoordB(:,2));
    set(gca,'fontsize',18); box on; axis equal; axis tight; view(2); set(gca,'YDir','reverse');
    title('Tracked displacements','fontweight','normal');
    xlabel(''); ylabel(''); cb = colorbar; set(cb,'fontsize',18);

    
    disp('Press any key to keep running code!')
    pause;

    %%%%% Store results %%%%%
    uv_B2A_prev{ImgSeqNum-1} = u_B2A_curr_refB;  % incremental displacement
    track_A2B_prev{ImgSeqNum-1} = track_A2B;
    track_B2A_prev{ImgSeqNum-1} = track_B2A;

end





