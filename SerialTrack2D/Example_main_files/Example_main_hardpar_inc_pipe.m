%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SerialTrack execute main file
% ===================================================
% Dimension:            2D
% Particle rigidity:    hard 
% Tracking mode:        incremental
% Syn or Exp:           exp
% Deformation mode:     pipe flow
%
% ===================================================
% Author: Jin Yang, Ph.D.
% Email: jyang526@wisc.edu -or-  aldicdvc@gmail.com 
% Date: 02/2022; 07/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% Initialization
close all; clear all; clc; clearvars -global
disp('************************************************');
disp('*** Welcome to SerialTrack Particle Tracking ***');
disp('************************************************');
addpath( './function/','./src/'); 

 
%% user defined parameters %%%%%

%%%%% Problem dimension and units %%%%%
MPTPara.DIM = 2;                % problem dimension
MPTPara.xstep = 0.075;          % unit: mm/px
MPTPara.tstep = 7.6923e-4;      % unit: s

%%%%% Code mode %%%%%
MPTPara.mode = 'inc'; % {'inc': incremental mode; 
                      %  'accum': cumulative mode; 
                      %  'dbf': double frame}

%%%%% Particle rigidity %%%%%
MPTPara.parType = 'hard'; % {'hard': hard particle; 
                          %  'soft': soft particle}
 
disp('************************************************');
disp(['Dimention: ',num2str(MPTPara.DIM)]);
disp(['Tracking mode: ',MPTPara.mode]);
disp(['Particle type: ',MPTPara.parType]);
disp('************************************************'); fprintf('\n');

%%%%% Image binary mask file %%%%%
MaskFileLoadingMode = 3; % {0: No mask file
                         %  1: Load only one mask file for all frames; 
                         %  2: Load one mask file for each single frame;
                         %  3: Load a MATLAB mat file for all frames;

if MaskFileLoadingMode == 3
    im_roi_mask_file_path = '.\img_par2track_pipe\im_roi.mat'; % TODO: Path of the mat file to be used as the mask file
else
    im_roi_mask_file_path = '';  
end
 

%%%%% Particle detection and localization parameters %%%%%
%%%%% Bead Parameters %%%%%
BeadPara.thres = 0.5;           % Threshold for detecting particles
BeadPara.beadRad = 3;           % Estimated radius of a single particle [px]
BeadPara.minSize = 2;           % Minimum area of a single particle [px^2]
BeadPara.maxSize = 20;          % Maximum area of a single particle [px^2]
BeadPara.winSize = [5, 5];      % Default [not used in 2D]
BeadPara.dccd = [1,1];          % Default [not used in 2D]
BeadPara.abc = [1,1];           % Default [not used in 2D]
BeadPara.forloop = 1;           % Default [not used in 2D]
BeadPara.randNoise = 1e-7;      % Default [not used in 2D]
BeadPara.PSF = [];              % PSF function; Example: PSF = fspecial('disk', BeadPara.beadRad-1 ); % Disk blur
BeadPara.color = 'white';       % Bead color: 'white' -or- 'black'


%% SerialTrack particle tracking

%%%%% Multiple particle tracking (MPT) Parameter %%%%%
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
MPTPara.distMissing = 2;       % Distance threshold to check whether particle has a match or not 


%%%% Postprocessing: merge trajectory segments %%%%%
distThres = 1;              % distance threshold to connect split trajectory segments [px]
extrapMethod = 'pchip';     % extrapolation scheme to connect split trajectory segments
                            % suggestion: 'nearest' for Brownian motion                          
minTrajSegLength = 10;      % the minimum length of trajectory segment that will be extrapolated [px]
maxGapTrajSeqLength = 0;    % the max frame# gap between connected trajectory segments


%%%%% Execute SerialTrack particle tracking %%%%%
if strcmp(MPTPara.mode,'inc')==1
    run_Serial_MPT_2D_hardpar_inc;
elseif strcmp(MPTPara.mode,'accum')==1
    run_Serial_MPT_2D_hardpar_accum;    
elseif strcmp(MPTPara.mode,'dbf')==1
    run_Serial_MPT_2D_hardpar_dbf;
end




