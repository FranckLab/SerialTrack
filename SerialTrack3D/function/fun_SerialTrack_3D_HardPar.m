function [parCoordB,uvw_B2A_refB,resultDisp,resultDefGrad,track_A2B,track_B2A] = fun_SerialTrack_3D_HardPar(...
    ImgSeqNum,ImgDef,BeadPara,MPTPara,parCoordA,parCoordB_prev,uvw_B2A_refB_prev)
%FUNCTION fun_SerialTrack_3D_HardPar
% Objective: to track 3D hard particle pairs 
% -----------------------------------------------
%
%       Input variables     Physical meaning
%       ================    ===================
%       ImgSeqNum           Image frame index
%       ImgDef              Current deformed image
%       BeadPara            Particle (bead) detection parameters
%       MPTPara             Particle linking parameters
%       parCoordA           Particle coordinates in the first frame
%       parCoordB_prev      Stored particle coordinates in previous frames
%       uvw_B2A_refB_prev   Stored particle centroid displacements in prev. frames
%
%       Output variables    Physical meaning
%       ================    ===================
%       parCoordB           Detected particles in current, deformed image
%       uvw_B2A_refB        Tracked displacement field based on current configuration
%       resultDisp          Tracked displacement vector 
%       resultDefGrad       Tracked displacement gradient vector
%       track_A2B           Tracked particle index links: parCoordA -> parCoordB
%       track_B2A           Tracked particle index links: parCoordB -> parCoordA
%
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

warning('off');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 1: Particle detection and localization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Load current deformed frame %%%%%
currImg = ImgDef(MPTPara.gridxyzROIRange.gridx(1):MPTPara.gridxyzROIRange.gridx(2), ...
                 MPTPara.gridxyzROIRange.gridy(1):MPTPara.gridxyzROIRange.gridy(2), ...
                 MPTPara.gridxyzROIRange.gridz(1):MPTPara.gridxyzROIRange.gridz(2));

%%%%% If PSF is non-empty, perform deconvolution %%%%%
if ~isempty(BeadPara.PSF)
    currImg = deconvlucy(currImg,BeadPara.PSF,6);
    disp('----- Deconvolution frame #',num2str(ImgSeqNum),' ------');
end

%%%%% Pre-process bead image if bead color is "black" %%%%%
if strcmp(BeadPara.color,'black')
    disp('Not implement yet.')
else
    currImg2 = currImg;
end

%%%%% Several methods to detect and localize particles %%%%%
try 
    BeadPara.detectionMethod = BeadPara.detectionMethod;
catch
    BeadPara.detectionMethod = 2;
end
% ----------------------------
%%%%% Method 1: TPT code %%%%%
if BeadPara.detectionMethod == 1 
    x{1}{ImgSeqNum} = locateParticles(double(currImg2)/max(double(currImg2(:))),BeadPara); % Detect particles
    x{1}{ImgSeqNum} = radialcenter3dvec(double(currImg2),x{1}{ImgSeqNum},BeadPara); % Localize particles
% ----------------------------
%%%%% Method 2: Modified TracTrac code %%%%%
elseif BeadPara.detectionMethod == 2
    x{1}{ImgSeqNum} = f_detect_particles3(double(currImg2)/max(double(currImg2(:))),BeadPara); %detect and localize particles
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Add MPTPara.gridxyzROIRange left-bottom corner point coordinates
x{1}{ImgSeqNum} = x{1}{ImgSeqNum} + [MPTPara.gridxyzROIRange.gridx(1)-1, ...
                                     MPTPara.gridxyzROIRange.gridy(1)-1, ...
                                     MPTPara.gridxyzROIRange.gridz(1)-1];
parCoordB = x{1}{ImgSeqNum};

% Remove bad coordinates that are out of image ROI
for tempi=1:3, parCoordB( parCoordB(:,tempi)>size(ImgDef,tempi), : ) = []; end
for tempi=1:3, parCoordB( parCoordB(:,tempi)<1, : ) = []; end


%%%%% Plot detected particles %%%%%
% close all;
% figure, plot3(parCoordA(:,1),parCoordA(:,2),parCoordA(:,3),'ro'); 
% view(2); box on; axis equal; axis tight; set(gca,'fontsize',18); 
% title('Detected particles in ref image','fontweight','normal');
% figure, plot3(parCoordB(:,1),parCoordB(:,2),parCoordB(:,3),'ro'); 
% view(2); box on; axis equal; axis tight; set(gca,'fontsize',18); 
% title('Detected particles in def image','fontweight','normal');


%%%%% Report detected beads # %%%%%
% disp(['Detected particle # in ref image: ',num2str(size(parCoordA,1))]);
disp(['Detected particle # in def image: ',num2str(size(parCoordB,1))]);
% disp('%%%%%% Detect particles: Done! %%%%%%'); fprintf('\n');



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
  
[matches_A2B,u_B2A_curr_refB,track_A2B] = f_track_serial_match3D( parCoordA, parCoordB, ...
   'f_o_s',MPTPara.f_o_s, 'n_neighborsMax', MPTPara.n_neighborsMax, 'n_neighborsMin', MPTPara.n_neighborsMin, ...
   'gbSolver', MPTPara.gbSolver, 'smoothness', MPTPara.smoothness, ...
   'outlrThres', MPTPara.outlrThres, 'maxIterNum', MPTPara.maxIterNum, ...
   'iterStopThres', MPTPara.iterStopThres, 'usePrevResults', MPTPara.usePrevResults, ...
   'strain_n_neighbors',MPTPara.strain_n_neighbors, 'strain_f_o_s',MPTPara.strain_f_o_s, ...
   'gridxyzROIRange',MPTPara.gridxyzROIRange, 'parCoordB_prev',parCoordB_prev, ...
   'uvw_B2A_prev',uvw_B2A_refB_prev, 'ImgSeqNum',ImgSeqNum, ...
   'MPTParaDistMissing',MPTPara.distMissing);


%%%%% Compute track_B2A %%%%%
matches_A2B = matches_A2B(matches_A2B(:,2)>0,:); % Build matches_A2B_temp
track_B2A = zeros(size(parCoordB,1), 1);
for tempi = 1:size(matches_A2B)
    track_B2A(matches_A2B(tempi,2)) = matches_A2B(tempi,1);
end  

 
%% %%%%% Plotting %%%%%%
% Compute displacement from tracked particles on deformed frame
disp_A2B_parCoordB = -u_B2A_curr_refB;
figure, plotCone3(parCoordB(:,1),parCoordB(:,2),parCoordB(:,3), ...
        disp_A2B_parCoordB(:,1),disp_A2B_parCoordB(:,2),disp_A2B_parCoordB(:,3));
set(gca,'fontsize',18); box on; axis equal; axis tight; view(3);  
title('Tracked displacements','fontweight','normal');
xlabel(''); ylabel(''); cb = colorbar; set(cb,'fontsize',18);



%% %%%%% Compute F = def grad = grad_u = (grad_x X - I) %%%%%

% strain_f_o_s: size of virtual strain gauge
% strain_n_neighbors: # of neighboring particles used in strain gauge

%%%%% Strain based on Moving Least Square Fitting in deformed configuration %%%%%
[XYZ_refB,U_B2A_refB,F_B2A_refB] = funCompDefGrad3(disp_A2B_parCoordB, parCoordB, MPTPara.strain_f_o_s, MPTPara.strain_n_neighbors);

%%%%% Strain based on Moving Least Square Fitting in reference configuration %%%%%
[XYZ_refA,U_A2B_refA,F_A2B_refA] = funCompDefGrad3(disp_A2B_parCoordB, parCoordB-disp_A2B_parCoordB, MPTPara.f_o_s, MPTPara.strain_n_neighbors);

 

%% %%%%% Store results %%%%%
uvw_B2A_refB = u_B2A_curr_refB;

resultDisp = struct('parCoordA',parCoordA,'parCoordB',parCoordB,'track_A2B', ...
                         track_A2B,'disp_A2B_parCoordB',disp_A2B_parCoordB);
                     
resultDefGrad = struct('XY_refA',XYZ_refA,'U_A2B_refA',U_A2B_refA,'F_A2B_refA',F_A2B_refA, ...
                        'XY_refB',XYZ_refB,'U_B2A_refB',U_B2A_refB,'F_B2A_refB',F_B2A_refB);






