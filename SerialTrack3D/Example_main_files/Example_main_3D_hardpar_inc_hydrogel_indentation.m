%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SerialTrack execute main file
% ===================================================
% Dimension:            3D
% Particle rigidity:    hard 
% Tracking mode:        incremental
% Syn or Exp:           exp
% Deformation mode:     hydrogel indentation
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

 
%% user defined parameters %%%%%

%%%%% Problem dimension and units %%%%%
MPTPara.DIM = 3;      % problem dimension
MPTPara.xstep = 0.4;  % unit: um/px
MPTPara.ystep = 0.4;  % unit: um/px
MPTPara.zstep = 0.4;  % unit: um/px
MPTPara.tstep = 1;    % unit: us

%%%%% Code mode %%%%%
MPTPara.mode = 'inc'; % {'inc': incremental mode; 
                      %  'accum': accumulative mode}

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
fileNameAll = 'vol_hydrogel_ind_20190504_cut_0*.mat';
fileFolder = './imgFolder/img_hydrogel_indentation/';

%%%%% Image binary mask file %%%%%
im_roi_mask_file_path = ''; % TODO: leave it as empty if there is no mask file


%%%%% Particle detection and localization parameters %%%%%

%%%%% particle detection and localization method %%%%%
BeadPara.detectionMethod = 2;   % Particle detection method: 1 = TPT (blob finding + radial projection), 
%                                                            2 = TracTrac (LoG blob finding + lsq fit of gaussian)

%%%%% Bead Parameter %%%%%
BeadPara.thres = 0.1;           % Threshold for detecting particles
BeadPara.beadRad = 3;           % Estimated radius of a single particle
BeadPara.minSize = 3;           % Minimum volume of a single particle [px^3]
BeadPara.maxSize = 20;          % Maximum volume of a single particle [px^3]
BeadPara.winSize = [5,5,5];     % Default [not used for method 2]
BeadPara.dccd = [1,1,1];        % Default [not used for method 2]
BeadPara.abc = [1,1,1];         % Default [not used for method 2]
BeadPara.forloop = 1;           % Default [not used for method 2]
BeadPara.randNoise = 1e-7;      % Default small number for background noise
BeadPara.PSF = [];              % PSF function; Example: PSF = fspecial('disk', BeadPara.beadRad-1 ); % Disk blur
BeadPara.color = 'white';       % Foreground (particle) color: options, 'white' or 'black'


%% SerialTrack particle tracking

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


%%%%% Execute SerialTrack particle tracking %%%%%
if strcmp(MPTPara.mode,'inc')==1
    run_Serial_MPT_3D_hardpar_inc;
elseif strcmp(MPTPara.mode,'accum')==1
    run_Serial_MPT_3D_hardpar_accum;    
end




