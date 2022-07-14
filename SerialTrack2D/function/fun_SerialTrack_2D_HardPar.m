function [parCoordB,uv_B2A_refB,resultDisp,resultDefGrad,track_A2B,track_B2A] = fun_SerialTrack_2D_HardPar(...
    ImgSeqNum,defImg,BeadPara,MPTPara,parCoordA,parCoordB_prev,uv_B2A_refB_prev)
%FUNCTION fun_SerialTrack_2D_HardPar
% Objective: to track 2D hard particle pairs 
% -----------------------------------------------
%
%       Input variables     Physical meaning
%       ================    ===================
%       ImgSeqNum           Image frame index
%       defImg              Current deformed image
%       BeadPara            Particle (bead) detection parameters
%       MPTPara             Particle linking parameters
%       parCoordA           Particle coordinates in the first frame
%       parCoordB_prev      Stored particle coordinates in previous frames
%       uv_B2A_refB_prev    Stored particle centroid displacements in prev. frames
%
%       Output variables    Physical meaning
%       ================    ===================
%       parCoordB           Detected particles in current, deformed image
%       uv_B2A_refB         Tracked displacement field based on current configuration
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
% SECTION 1: Particle detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Load current deformed frame %%%%%
currImg = defImg(MPTPara.gridxyROIRange.gridx(1):MPTPara.gridxyROIRange.gridx(2), ...
    MPTPara.gridxyROIRange.gridy(1):MPTPara.gridxyROIRange.gridy(2));

%%%%% If PSF is non-empty, perform deconvolution %%%%%
if ~isempty(BeadPara.PSF)
    currImg = deconvlucy(currImg,BeadPara.PSF,6);
    disp('----- Deconvolution frame #',num2str(ImgSeqNum),' ------');
end


%%%%% Pre-process bead image if bead color is "black" %%%%%
if strcmp(BeadPara.color,'black')
%     ImgGauss = imgaussfilt(imgaussfilt(currImg,1),1);
%     ImgGauss(ImgGauss > BeadPara.thres*max(double(currImg(:)))) = 0;
%     bw = imbinarize(uint8(ImgGauss),'adaptive','ForegroundPolarity','dark','Sensitivity',0.8); % figure, imshow(bws2);
%     bws2 = bwareaopen(bw,round(pi*BeadPara.minSize^2)); % remove all object containing fewer than BeadPara.minSize
%     removeobjradius = BeadPara.minSize; % fill a gap in the pen's cap
%     se = strel('disk',removeobjradius);
%     bws2 = imclose(bws2,se);
    currImg_norm = double(currImg)/max(double(currImg(:)));
    currImg2_norm = imcomplement(currImg_norm); % figure, imshow(uint8(currImg2));
else
    currImg2_norm = double(currImg)/max(double(currImg(:)));
end


%%%%% Detect particles %%%%%
x{1}{ImgSeqNum} = f_detect_particles(currImg2_norm,BeadPara);

% Add MPTPara.gridxyROIRange left-bottom corner point coordinates
x{1}{ImgSeqNum} = x{1}{ImgSeqNum} + [MPTPara.gridxyROIRange.gridx(1)-1, MPTPara.gridxyROIRange.gridy(1)-1];
parCoordB = x{1}{ImgSeqNum};

% Remove bad coordinates that are out of image ROI
for tempi=1:2, parCoordB( parCoordB(:,tempi)>size(defImg,tempi), : ) = []; end
for tempi=1:2, parCoordB( parCoordB(:,tempi)<1, : ) = []; end


%%%%% Plot detected particles %%%%%
% close all; figure,
% hold on; plot( parCoordA(:,1), parCoordA(:,2), 'bo');
% view(2); box on; axis equal; axis tight; set(gca,'fontsize',18);
% title('Detected particles in ref image','fontweight','normal');
% 
% % figure,
% hold on; plot(parCoordB(:,1), parCoordB(:,2),'ro');
% view(2); box on; axis equal; axis tight; set(gca,'fontsize',18);
% title('Detected particles in defor image','fontweight','normal');
% pause;


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
  
[matches_A2B,u_B2A_curr_refB,track_A2B] = f_track_serial_match2D( parCoordA, parCoordB, ...
   'f_o_s',MPTPara.f_o_s, 'n_neighborsMax', MPTPara.n_neighborsMax, 'n_neighborsMin', MPTPara.n_neighborsMin, ...
   'locSolver',MPTPara.locSolver,'gbSolver', MPTPara.gbSolver, 'smoothness', MPTPara.smoothness, ...
   'outlrThres', MPTPara.outlrThres, 'maxIterNum', MPTPara.maxIterNum, ...
   'iterStopThres', MPTPara.iterStopThres, 'usePrevResults', MPTPara.usePrevResults, ...
   'strain_n_neighbors',MPTPara.strain_n_neighbors, 'strain_f_o_s',MPTPara.strain_f_o_s, ...
   'gridxyROIRange',MPTPara.gridxyROIRange, 'parCoordB_prev',parCoordB_prev, ...
   'uv_B2A_prev',uv_B2A_refB_prev, 'ImgSeqNum',ImgSeqNum, ...
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
figure, plotCone2(parCoordB(:,1),parCoordB(:,2),disp_A2B_parCoordB(:,1),disp_A2B_parCoordB(:,2));
set(gca,'fontsize',18); box on; axis equal; axis tight; view(2); set(gca,'YDir','reverse');
title('Tracked displacements','fontweight','normal');
xlabel(''); ylabel(''); cb = colorbar; set(cb,'fontsize',18);



%% %%%%% Compute F = def grad = grad_u = (grad_x X - I) %%%%%

% strain_f_o_s: size of virtual strain gauge
% strain_n_neighbors: # of neighboring particles used in strain gauge

%%%%% Strain based on Moving Least Square Fitting in deformed configuration %%%%%
[XY_refB,U_B2A_refB,F_B2A_refB] = funCompDefGrad2(disp_A2B_parCoordB, parCoordB, MPTPara.strain_f_o_s, MPTPara.strain_n_neighbors);

%%%%% Strain based on Moving Least Square Fitting in reference configuration %%%%%
[XY_refA,U_A2B_refA,F_A2B_refA] = funCompDefGrad2(disp_A2B_parCoordB, parCoordB-disp_A2B_parCoordB, MPTPara.f_o_s, MPTPara.strain_n_neighbors);

 

%% %%%%% Store results %%%%%
uv_B2A_refB = u_B2A_curr_refB;

resultDisp = struct('parCoordA',parCoordA,'parCoordB',parCoordB,'track_A2B', ...
                         track_A2B,'disp_A2B_parCoordB',disp_A2B_parCoordB);
                     
resultDefGrad = struct('XY_refA',XY_refA,'U_A2B_refA',U_A2B_refA,'F_A2B_refA',F_A2B_refA, ...
                        'XY_refB',XY_refB,'U_B2A_refB',U_B2A_refB,'F_B2A_refB',F_B2A_refB);






