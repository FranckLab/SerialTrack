%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SerialTrack execute main file
% ===================================================
% Dimension:            3D
% Particle rigidity:    hard 
% Tracking mode:        incremental
% Syn or Exp:           syn
% Deformation mode:     rigid body motions including translations & rotations
%
% ===================================================
% Author: Jin Yang, Ph.D.
% Email: jyang526@wisc.edu -or-  aldicdvc@gmail.com 
% Date: 02/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization 
close all; clear all; clc; clearvars -global
disp('************************************************');
disp('*** Welcome to SerialTrack Particle Tracking ***');
disp('************************************************');
addpath( './function/','./src/','./Scatter2Grid3D' ); 

 
%% User defined parameters %%%%% 

%%%%% Problem dimension and units %%%%%
MPTPara.DIM = 3;    % problem dimension
MPTPara.xstep = 1;  % unit: um/px
MPTPara.tstep = 1;  % unit: us

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

%%%%% SerialTrack path %%%%%
SerialTrackPath = pwd; % % TODO: modify the path of "SerialTrack3D";
% Example: SerialTrackPath = 'D:\MATLAB\SerialTrack-main\SerialTrack3D';

%%%%% Volumetric image path %%%%%
% fileNameAll = 'VM_3D*.mat';
% fileFolder = '';

%%%%% Synthetic cases %%%%%
DefType = 'stretch';        % Loading type: {'translation','stretch','simpleshear','rotation'}
SeedingDensityType = 2;     % Particle seeding density: {1,2,3,4} --> {10,100,300,1000}*1e-6 particles per voxel 
fileNameAll = 'vol_*.mat';  % file name(s)

% ----- file folder name -----
if strcmp(DefType,'translation')==1
    fileFolder = ['./imgFolder/img_syn_hardpar/img_trans_hardpar_sd',num2str(SeedingDensityType)];
elseif strcmp(DefType,'rotation')==1
    fileFolder = ['./imgFolder/img_syn_hardpar/img_rot_hardpar_sd',num2str(SeedingDensityType)];
elseif strcmp(DefType,'stretch')==1
    fileFolder = ['./imgFolder/img_syn_hardpar/img_str_hardpar_sd',num2str(SeedingDensityType)];
elseif strcmp(DefType,'simpleshear')==1
    fileFolder = ['./imgFolder/img_syn_hardpar/img_simshear_hardpar_sd',num2str(SeedingDensityType)];
else
end

%%%%% Image binary mask file %%%%%
im_roi_mask_file_path = ''; % TODO: leave it as empty if there is no mask file


%%%%% Particle detection and localization parameters %%%%%

%%%%% Bead detection and localization method %%%%%
BeadPara.detectionMethod = 1;   % Particle detection method: 1 = TPT (blob finding + radial projection), 
%                                                            2 = TracTrac (LoG blob finding + lsq fit of gaussian)
%%%%% Bead Parameters %%%%%
BeadPara.thres = 0.5;           % Threshold for detecting particles
BeadPara.beadRad = 3;           % Estimated radius of a single particle [px]
BeadPara.minSize = 4;           % Minimum volume of a single particle [px^3]
BeadPara.maxSize = 100;         % Maximum volume of a single particle [px^3]
BeadPara.winSize = [5,5,5];     % Default window size for particle localization [used for method 1]
BeadPara.dccd = [1,1,1];        % Default [grid spacing for localization in method 1]
BeadPara.abc = [1,1,1];         % Default [grid spacing factor for localization in method 1]
BeadPara.forloop = 1;           % Default ["for" of linear indexing in method 1]
BeadPara.randNoise = 1e-7;      % Default [small amount of background noise method 1]
BeadPara.PSF = [];              % PSF function; Example: PSF = fspecial('disk', BeadPara.beadRad-1 ); % Disk blur
BeadPara.color = 'white';       % Foreground (particle) color: options, 'white' or 'black'


%% SerialTrack particle tracking

%%%%% Multiple particle tracking (MPT) Parameters %%%%%
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
MPTPara.distMissing = 5;         % Distance threshold to check whether particle has a match or not [px]


%%%% Postprocessing: merge trajectory segments %%%%%
distThres = 1;              % distance threshold to connect split trajectory segments
extrapMethod = 'pchip';     % extrapolation scheme to connect split trajectory segments
                            % suggestion: 'nearest' for Brownian motion                          
minTrajSegLength = 10;      % the minimum length of trajectory segment that will be extrapolate 
maxGapTrajSeqLength = 0;    % the max frame# gap between connected trajectory segments

 
%%%%% Execute SerialTrack particle tracking %%%%%
if strcmp(MPTPara.mode,'inc')==1
    run_Serial_MPT_3D_hardpar_inc;
elseif strcmp(MPTPara.mode,'accum')==1
    run_Serial_MPT_3D_hardpar_accum;    
end
 

