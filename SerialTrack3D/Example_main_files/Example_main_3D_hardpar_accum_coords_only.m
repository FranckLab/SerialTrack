%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SerialTrack execute main file
% ===================================================
% Dimension:            3D
% Particle rigidity:    hard
% Tracking mode:        cumulative
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
addpath( './function/','./src/','./Scatter2Grid3D' );

%%%%% SerialTrack path %%%%%
SerialTrackPath = pwd; % % TODO: modify the path of "SerialTrack3D";
% Example: SerialTrackPath = 'D:\MATLAB\SerialTrack-main\SerialTrack3D';


%% Load tracked particle coordinates "parCoordA" and "parCoordB"

fprintf('\n'); disp('************************************************');
disp('This example is to demo how SerialTrack links');
disp('particles whose centroids are already provided.');
disp('************************************************'); fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Option I: If you only have two images %%%%%
% % TODO: load particle centroid locations in the ref image as "parCoordA"
% load('.\Example_main_files\data_hydrogel_ind_parCoordA.mat');
% % TODO: load particle centroid locations in the def image as "parCoordB"
% load('.\Example_main_files\data_hydrogel_ind_parCoordB.mat');
% 
% parCoord_prev{1} = parCoordA;
% parCoord_prev{2} = parCoordB;

%%%%% Option II: If you have more than two images %%%%%
% TODO: load particle centroid locations and name it as a cell array "parCoord_prev"
load('.\Example_main_files\data_syn_simple_shear_parCoord.mat');
parCoordA = parCoord_prev{1}; 

disp('%%%%%% Load tracked particle corodinates: Done! %%%%%%'); fprintf('\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Image binary mask file %%%%%
im_roi_mask_file_path = ''; % TODO: leave it as empty if there is no mask file

%%%%% Load image mask file %%%%%
try load(im_roi_mask_file_path); catch; end
try MPTPara.ImgRefMask = im_roi'; % Load stored image roi if existed
catch, MPTPara.ImgRefMask = ones(max(ceil(parCoordA(:,:)))); % Set up default image mask file
end
disp('%%%%%% Load image mask file: Done! %%%%%%'); fprintf('\n');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 1: User defined parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Update MPTPara %%%%%
MPTPara.gridxyzROIRange.gridx = [min(floor(parCoordA(:,1))), max(ceil(parCoordA(:,1)))];
MPTPara.gridxyzROIRange.gridy = [min(floor(parCoordA(:,2))), max(ceil(parCoordA(:,2)))];
MPTPara.gridxyzROIRange.gridz = [min(floor(parCoordA(:,3))), max(ceil(parCoordA(:,3)))];

%%%%% Problem dimension and units %%%%%
MPTPara.DIM = 3;      % problem dimension
MPTPara.xstep = 0.4;  % unit: um/px
MPTPara.ystep = 0.4;  % unit: um/px
MPTPara.zstep = 0.4;  % unit: um/px
MPTPara.tstep = 1;    % unit: us

%%%%% Code mode %%%%%
MPTPara.mode = 'accum'; % {'inc': incremental mode;
                        %  'accum': cumulative mode}

%%%%% Particle rigidity %%%%%
MPTPara.parType = 'hard'; % {'hard': hard particle;
                          %  'soft': soft particle}

disp('************************************************');
disp(['Dimention: ',num2str(MPTPara.DIM)]);
disp(['Tracking mode: ',MPTPara.mode]);
disp(['Particle type: ',MPTPara.parType]);
disp('************************************************'); fprintf('\n');

%%%%% Multiple particle tracking (MPT) Parameter %%%%%
MPTPara.f_o_s = 60;              % Size of search field: max(|u|,|v|,|w|) [px]
MPTPara.n_neighborsMax = 25;     % Max # of neighboring particles
MPTPara.n_neighborsMin = 1;      % Min # of neighboring particles
MPTPara.gbSolver = 2;            % Global step solver: 1-moving least square fitting; 2-global regularization; 3-ADMM iterations
MPTPara.smoothness = 1e-1;       % Coefficient of regularization
MPTPara.outlrThres = 5;          % Threshold for removing outliers in MPT
MPTPara.maxIterNum = 20;         % Max ADMM iteration number
MPTPara.iterStopThres = 1e-3;    % ADMM iteration stopping threshold
MPTPara.strain_n_neighbors = 20; % # of neighboring particles used in strain gauge
MPTPara.strain_f_o_s = 60;       % Size of virtual strain gauge [px]
MPTPara.usePrevResults = 0;      % Whether use previous results or not: 0-no; 1-yes;
MPTPara.distMissing = 5;         % Distance threshold to check whether particle has a match or not


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
uvw_B2A_prev = cell(length(parCoord_prev)-1,1);

for ImgSeqNum = 2 : length(parCoord_prev) % "ImgSeqNum" is the frame index

    disp(['====== Frame #',num2str(ImgSeqNum),' ======']);

    [matches_A2B,uvw_B2A_curr_refB,track_A2B] = f_track_serial_match3D( ...
        parCoord_prev{1}, parCoord_prev{ImgSeqNum}, ...
        'f_o_s',MPTPara.f_o_s, 'n_neighborsMax', MPTPara.n_neighborsMax, 'n_neighborsMin', MPTPara.n_neighborsMin, ...
        'gbSolver', MPTPara.gbSolver, 'smoothness', MPTPara.smoothness, ...
        'outlrThres', MPTPara.outlrThres, 'maxIterNum', MPTPara.maxIterNum, ...
        'iterStopThres', MPTPara.iterStopThres, 'usePrevResults', MPTPara.usePrevResults, ...
        'strain_n_neighbors',MPTPara.strain_n_neighbors, 'strain_f_o_s',MPTPara.strain_f_o_s, ...
        'gridxyzROIRange',MPTPara.gridxyzROIRange, 'parCoordB_prev',parCoord_prev, ...
        'uvw_B2A_prev',uvw_B2A_prev, 'ImgSeqNum',ImgSeqNum, ...
        'MPTParaDistMissing',MPTPara.distMissing);


    %%%%% Compute track_B2A %%%%%
    matches_A2B = matches_A2B(matches_A2B(:,2)>0,:); % Build matches_A2B_temp
    track_B2A = zeros(size(parCoord_prev{ImgSeqNum},1), 1);
    for tempi = 1:size(matches_A2B)
        track_B2A(matches_A2B(tempi,2)) = matches_A2B(tempi,1);
    end


    %% %%%%% Plotting %%%%%%
    % Compute displacement from tracked particles on deformed frame
    disp_A2B_parCoordB = -uvw_B2A_curr_refB;
    figure, plotCone3(parCoord_prev{ImgSeqNum}(:,1),parCoord_prev{ImgSeqNum}(:,2),parCoord_prev{ImgSeqNum}(:,3), ...
        disp_A2B_parCoordB(:,1),disp_A2B_parCoordB(:,2),disp_A2B_parCoordB(:,3));
    set(gca,'fontsize',18); box on; axis equal; axis tight; view(3);
    title('Tracked displacements','fontweight','normal');
    xlabel(''); ylabel(''); cb = colorbar; set(cb,'fontsize',18);


    %%%%% Store results %%%%%
    uvw_B2A_prev{ImgSeqNum-1} = uvw_B2A_curr_refB;  % incremental displacement
    track_A2B_prev{ImgSeqNum-1} = track_A2B;
    track_B2A_prev{ImgSeqNum-1} = track_B2A;

end





