function [xyGrid_prevCurr,uvGrid_B2A_refB_prevCurr,track_A2B,parCoordBCurr] = fun_SerialTrack_2D_SoftPar(...
    ImgSeqNum,defImg,BeadPara,MPTPara,parCoordA,xyGrid_prev,uvGrid_B2A_refB_prev)
%FUNCTION fun_SerialTrack_2D_SoftPar
% Objective: to track 2D soft particle pairs
% -----------------------------------------------
%
%       Input variables         Physical meaning
%       ================        ===================
%       ImgSeqNum               Image frame index
%       defImg                  Current deformed image
%       BeadPara                Particle (bead) detection parameters
%       MPTPara                 Particle linking parameters
%       parCoordA               Particle coordinates in the first frame
%       xyGrid_prev             Stored point coordinates in previous frames
%       uvGrid_B2A_refB_prev    Stored displacement results in prev. frames
%
%       Output variables    Physical meaning
%       ================    ===================
%       xyGrid_prevCurr           xy-grid nodal points coordinates in current frame
%       uvGrid_B2A_refB_prevCurr  Tracked displacement field based on current xy-grid mesh
%       track_A2B           Tracked particle index links: parCoordA -> parCoordB
%       parCoordB           Detected particles in current, deformed image
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



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Initialize iteration %%%%%
warning('off'); iterNum = 0;

try f_o_s = MPTPara.f_o_s;                   catch f_o_s = 20;              end   % Size of search field [px]
try n_neighborsMax = MPTPara.n_neighborsMax; catch n_neighborsMax = 25;     end   % Max # of neighboring particles
try n_neighborsMin = MPTPara.n_neighborsMin; catch n_neighborsMin = 1;      end   % Min # of neighboring particles
try locSolver = MPTPara.locSolver;           catch locSolver = 1;           end   % Local solver: 1-topology-based feature; 2-histogram-based feature first and then topology-based feature;
try gbSolver = MPTPara.gbSolver;             catch gbSolver = 1;            end   % Global step solver: 1-moving least square fitting; 2-global regularization; 3-ADMM iterations
try smoothness = MPTPara.smoothness;         catch smoothness = 1e-2;       end   % Coefficient of regularization
try outlrThres = MPTPara.outlrThres;         catch outlrThres = 2;          end   % Threshold for removing outliers in TPT (Westerweel et al., Exp.Fluids, 2005)
try maxIterNum = MPTPara.maxIterNum;         catch maxIterNum = 20;         end   % Max ADMM iteration number
try iterStopThres = MPTPara.iterStopThres;   catch iterStopThres = 1e-2;    end   % ADMM iteration stopping threshold
try usePrevResults = MPTPara.usePrevResults; catch usePrevResults = 0;      end   % Whether use previous results or not
try strain_n_neighbors = MPTPara.strain_n_neighbors;
catch strain_n_neighbors = 20; end   % Size of virtual strain gauge [px]
try strain_f_o_s = MPTPara.strain_f_o_s;     catch strain_f_o_s = 50;       end   % # of neighboring particles used in strain gauge
try gridxyROIRange = MPTPara.gridxyROIRange;
catch
    gridxyROIRange.gridx = [min(parCoordA(:,1)),max(parCoordA(:,1))]; % ROI-gridx
    gridxyROIRange.gridy = [min(parCoordA(:,2)),max(parCoordA(:,2))]; % ROI-gridy
end

sxy = min([round(0.5*f_o_s),20])*[1,1];   % Grid size for regularization


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Load current deformed frame %%%%%
currImg = defImg(MPTPara.gridxyROIRange.gridx(1):MPTPara.gridxyROIRange.gridx(2), ...
    MPTPara.gridxyROIRange.gridy(1):MPTPara.gridxyROIRange.gridy(2));

%%%%% If PSF is non-empty, perform deconvolution %%%%%
if ~isempty(BeadPara.PSF)
    currImg = deconvlucy(currImg,BeadPara.PSF,6);
    disp('----- Deconvolution frame #',num2str(ImgSeqNum),' ------');
end

%%%%% Initialize two grid meshes: "fine-mesh is the image pixel mesh"
[xGrid_Img,yGrid_Img] = ndgrid( gridxyROIRange.gridx(1):gridxyROIRange.gridx(2), ...
    gridxyROIRange.gridy(1):gridxyROIRange.gridy(2) );
[xGrid,yGrid] = ndgrid(gridxyROIRange.gridx(1) : sxy(1) : gridxyROIRange.gridx(2), ...
    gridxyROIRange.gridy(1) : sxy(2) : gridxyROIRange.gridy(2));

uGrid_B2A_refB = 0*xGrid; vGrid_B2A_refB = 0*xGrid;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ADMM iteration while loop %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while iterNum < maxIterNum
    
    %%%%% # of neighboring particles: (topology-based search) --> 1 (nearest neighboring search) %%%%%
    n_neighbors = round( n_neighborsMin + exp(-0.5*iterNum)*(n_neighborsMax-n_neighborsMin) );
    iterNum = iterNum+1; disp(['------ Iter #',num2str(iterNum),' ------']); % current iteration #
    
    %%%%% Warp deformed image to re-detect beads %%%%%
    if iterNum == 1
        
        %%%%% Initialize "xGridCurr_Img_refB" and "yGridCurr_Img_refB" %%%%%
        xGridCurr_Img_refB = xGrid_Img;
        yGridCurr_Img_refB = yGrid_Img;
        
        %%%%% Use previous results to estimate the disp. in a new frame %%%%%
        if ImgSeqNum > 2 && usePrevResults
            
            %%%%% POD prediction %%%%%
            [tempu,tempv] = funInitGuess(xyGrid_prev,uvGrid_B2A_refB_prev,[xGrid(:),yGrid(:)],ImgSeqNum);
            
            uGrid_B2A_refB = reshape(tempu,size(xGrid));
            vGrid_B2A_refB = reshape(tempv,size(xGrid));
            
            Fx = scatteredInterpolant([xGrid(:),yGrid(:)],tempu,'linear','linear');
            tempu = Fx(xGrid_Img(:),yGrid_Img(:));
            Fy = scatteredInterpolant([xGrid(:),yGrid(:)],tempv,'linear','linear');
            tempv = Fy(xGrid_Img(:),yGrid_Img(:));
            
            %%%%% Add to current displacement vectors %%%%%
            xGridCurr_Img_refB = xGridCurr_Img_refB - reshape(tempu,size(xGrid_Img,1),size(xGrid_Img,2));
            yGridCurr_Img_refB = yGridCurr_Img_refB - reshape(tempv,size(xGrid_Img,1),size(xGrid_Img,2));
            
            %%%%% Warp deformed image %%%%%
            % --- MATLAB default interpolation function ---
            % currImg = interp2(ImgDef,yGridCurr_Img_refB,xGridCurr_Img_refB,'spline');
            % --- External "ba_interp2" code which runs faster ---
            currImg = ba_interp2(defImg,yGridCurr_Img_refB,xGridCurr_Img_refB,'cubic');
            
            %%%%% Remove outside grayscale values %%%%%
            [row1,~] = find( xGridCurr_Img_refB(:) < gridxyROIRange.gridx(1) );
            [row2,~] = find( xGridCurr_Img_refB(:) > gridxyROIRange.gridx(2) );
            [row3,~] = find( yGridCurr_Img_refB(:) < gridxyROIRange.gridy(1) );
            [row4,~] = find( yGridCurr_Img_refB(:) > gridxyROIRange.gridy(2) );
            row1234 = unique([row1;row2;row3;row4]);
            currImg(row1234) = nan;
            
        end %%%%% END of if
        
    elseif iterNum > 1
        
        %%%%% Iterative additive displacement %%%%%
        Fx = scatteredInterpolant([xGrid(:),yGrid(:)],uGrid_B2A_refB_iter(:),'linear','linear');
        tempu = Fx(xGrid_Img(:),yGrid_Img(:));
        Fy = scatteredInterpolant([xGrid(:),yGrid(:)],vGrid_B2A_refB_iter(:),'linear','linear');
        tempv = Fy(xGrid_Img(:),yGrid_Img(:));
        
        %%%%% Check convergence %%%%%
        %%%%% Norm of disp update %%%%%
        uv_B2A_updateNorm = sqrt( ( norm(tempu)^2+norm(tempv)^2 )/length(tempu) );
        disp(['Disp update norm: ',num2str(uv_B2A_updateNorm)]); % fprintf('\n');
        
        %%%%% Do at most five iterations that "matchRatio==1" %%%%%
        if matchRatio > 0.999
            matchRatioEqualsOneTime = matchRatioEqualsOneTime+1;
        end
        
        %%%%% Stopping criterion %%%%%
        if (uv_B2A_updateNorm<iterStopThres) || (matchRatioEqualsOneTime>5)
            disp(['----- Converged! ------']); fprintf('\n');
            break
        end
        
        %%%%% If not convergent yet,
        xGridCurr_Img_refB = xGridCurr_Img_refB - reshape(tempu,size(xGrid_Img,1),size(xGrid_Img,2));
        yGridCurr_Img_refB = yGridCurr_Img_refB - reshape(tempv,size(xGrid_Img,1),size(xGrid_Img,2));
        
        uGrid_B2A_refB = uGrid_B2A_refB + uGrid_B2A_refB_iter;
        vGrid_B2A_refB = vGrid_B2A_refB + vGrid_B2A_refB_iter;
        
        %%%%% Warp deformed image %%%%%
        currImg = ba_interp2(defImg,yGridCurr_Img_refB,xGridCurr_Img_refB,'cubic');
        
        %%%%% Remove outside grayscale values %%%%%
        [row1,~] = find( xGridCurr_Img_refB(:) < gridxyROIRange.gridx(1) );
        [row2,~] = find( xGridCurr_Img_refB(:) > gridxyROIRange.gridx(2) );
        [row3,~] = find( yGridCurr_Img_refB(:) < gridxyROIRange.gridy(1) );
        [row4,~] = find( yGridCurr_Img_refB(:) > gridxyROIRange.gridy(2) );
        
        row1234 = unique([row1;row2;row3;row4]);
        currImg(row1234) = nan;
        
    end
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SECTION 1: Particle detection and localization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    %%%%% Detect and localize particles %%%%%
    x{1}{ImgSeqNum} = f_detect_particles(currImg2_norm,BeadPara);
    
    %%%%% Add back gridxyROIRange %%%%%
    x{1}{ImgSeqNum} = x{1}{ImgSeqNum} + [gridxyROIRange.gridx(1)-1, gridxyROIRange.gridy(1)-1];
    parCoordBCurr = x{1}{ImgSeqNum};
    
    %%%%% Remove bad coordinates outside image ROI %%%%%
    for tempi=1:2, parCoordBCurr( parCoordBCurr(:,tempi)>size(defImg,tempi), : ) = []; end
    for tempi=1:2, parCoordBCurr( parCoordBCurr(:,tempi)<1, : ) = []; end
    
    %%%%% Plot %%%%%
    % figure, imshow(uint8(currImg)')
    % hold on; plot( parCoordBCurr(:,1)-(gridxyROIRange.gridx(1)-1), ...
    %     parCoordBCurr(:,2)-(gridxyROIRange.gridy(1)-1), 'r.');
    % view(2); box on; axis equal; axis tight; set(gca,'fontsize',18);
    % title('Detected particles in ref image','fontweight','normal');
    %
    % pause;
    
    %%%%% Report detected beads # %%%%%
    disp(['Detected particle # in def image: ',num2str(size(parCoordBCurr,1))]);
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SECTION 2: Remove non-paired particles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% Find particles that both exist in ref and def images %%%%%
    if iterNum == 1 || n_neighbors > 4 % n_neighborsMax/2  % Initialize ParNotMissing
        
        parNotMissingIndA = [1:size(parCoordA,1)]';
        parNotMissingIndBCurr = [1:size(parCoordBCurr,1)]'; % Remove missing particles in both ref and def images
        matchRatioEqualsOneTime = 0;
        
    else % iterNum > 1, find "ParNotMissing"
        
        %%%%% To detect missing particles in both reference and deformed images %%%%%
        neighborInd_BA = knnsearch(parCoordBCurr(:,1:2),parCoordA(:,1:2),'K',1); % Find pts in BCurr near pts in ACurr
        dist_BA = sqrt( sum((parCoordBCurr(neighborInd_BA,1:2) -  parCoordA).^2,2) ); % figure, h = histogram(dist_BCurrA);
        % [AParNotMissing,~] = find(dist_BCurrA < max([2, max(abs(u_B2A_curr(:)))])); % Find particles not missing in Particle A
        [parNotMissingIndA,~] = find(dist_BA < max([ MPTPara.distMissing, 0.5*BeadPara.beadRad])); % Find particles not missing in Particle A
        
        neighborInd_AB = knnsearch(parCoordA(:,1:2),parCoordBCurr(:,1:2),'K',1); % Find pts in ACurr near pts in BCurr
        dist_AB = sqrt( sum((parCoordA(neighborInd_AB,1:2) -  parCoordBCurr).^2,2) );
        % [BCurrParNotMissing,~] = find(dist_ABCurr < max([2, max(abs(u_B2A_curr(:)))])); % Find particles not missing in Particle BCurr
        [parNotMissingIndBCurr,~] = find(dist_AB < max([ MPTPara.distMissing, 0.5*BeadPara.beadRad]));  % Find particles not missing in Particle BCurr
        
    end
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SECTION 3: Particle tracking
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
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Local step: Neighbor topology match %%%%%
    matches_A2B = []; % Initialzation matches_A2B
    while isempty(matches_A2B) && n_neighborsMax < length(parNotMissingIndA)
        
        if locSolver == 1 && n_neighbors > 2 % topology-based feature
            matches_A2B = f_track_neightopo_match2( parCoordA(parNotMissingIndA,:), parCoordBCurr(parNotMissingIndBCurr,:), f_o_s, n_neighbors );
        elseif locSolver == 2 && n_neighbors > 2 % histogram-based feature
            matches_A2B = f_track_hist_match( parCoordA(parNotMissingIndA,:), parCoordBCurr(parNotMissingIndBCurr,:), f_o_s, n_neighbors, gauss_interp );
        else % nearest neighbor search
            matches_A2B = f_track_nearest_neighbour2( parCoordA(parNotMissingIndA,:), parCoordBCurr(parNotMissingIndBCurr,:), f_o_s );
        end
        if isempty(matches_A2B) == 1
            n_neighborsMax = round(n_neighborsMax + 5);
        else
            matches_A2B = [parNotMissingIndA(matches_A2B(:,1)), parNotMissingIndBCurr(matches_A2B(:,2))];
            [track_A2B, u_A2B] = funCompDisp2(parCoordA, parCoordBCurr, matches_A2B, outlrThres);
            matchRatio = size(matches_A2B,1) / length(parNotMissingIndA);
            disp( ['Tracking ratio: ', num2str( size(matches_A2B,1) ),'/', num2str(length(parNotMissingIndA)), ' = ', num2str(matchRatio) ]);
        end
        
    end
    if isempty(matches_A2B) == 1, disp('No matches are found. Maybe bad MPT parameters!'); pause; end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Global step: Compute kinematically compatible displacement %%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Method I: local moving least squares %%%%%
    if gbSolver==1
        
        [XY_B2A_refB,U_B2A_refB,F_B2A_refB] = funCompDefGrad2(-u_A2B, parCoordBCurr(track_A2B(track_A2B>0),:), f_o_s, n_neighbors);
        
        [row,~] = find(isnan(U_B2A_refB(1:2:end)) == 1);        % Find nans
        XY_B2A_refB(row,:) = [];                                % remove nans
        U_B2A_refB([2*row-1; 2*row ]) = [];                     % remove nans
        F_B2A_refB([4*row-3; 4*row-2; 4*row-1; 4*row]) = [];    % remove nans
        
        Fx = scatteredInterpolant(XY_B2A_refB,U_B2A_refB(1:2:end),'linear','linear');
        tempu = Fx(parCoordBCurr);
        Fy = scatteredInterpolant(XY_B2A_refB,U_B2A_refB(2:2:end),'linear','linear');
        tempv = Fy(parCoordBCurr);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Method II: global regularization %%%%%
    elseif gbSolver==2
        
        tempUV_B2A = -u_A2B;
        tempXY_refB = parCoordBCurr(track_A2B(track_A2B>0),:);
        
        try smoothness = smoothness; % Smoothness for regularization; "smoothness=0" means no regularization
        catch smoothness = 1e-3; % By default
        end
        try
            [~,~,uGrid_B2A_refB_iter] = funScatter2Grid2D(tempXY_refB(:,1),tempXY_refB(:,2),tempUV_B2A(:,1),sxy,smoothness,xGrid,yGrid);
            [~,~,vGrid_B2A_refB_iter] = funScatter2Grid2D(tempXY_refB(:,1),tempXY_refB(:,2),tempUV_B2A(:,2),sxy,smoothness,xGrid,yGrid);
        catch
            [~,~,uGrid_B2A_refB_iter] = funScatter2Grid2D(tempXY_refB(:,1),tempXY_refB(:,2),tempUV_B2A(:,1),sxy,0,xGrid,yGrid);
            [~,~,vGrid_B2A_refB_iter] = funScatter2Grid2D(tempXY_refB(:,1),tempXY_refB(:,2),tempUV_B2A(:,2),sxy,0,xGrid,yGrid);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Method III: augmented Lagrangian regularization %%%%%
    elseif gbSolver==3
        
        tempUV_B2A = -u_A2B;
        tempXY_refB = parCoordBCurr(track_A2B(track_A2B>0),:);
        try smoothness = smoothness; % Smoothness for regularization; "smoothness=0" means no regularization
        catch smoothness = 1e-3; % By default
        end
        try
            [~,~,uGrid_B2A_refB_iter] = funScatter2Grid2D(tempXY_refB(:,1),tempXY_refB(:,2),tempUV_B2A(:,1),sxy,smoothness,xGrid,yGrid);
            [~,~,vGrid_B2A_refB_iter] = funScatter2Grid2D(tempXY_refB(:,1),tempXY_refB(:,2),tempUV_B2A(:,2),sxy,smoothness,xGrid,yGrid);
        catch
            [~,~,uGrid_B2A_refB_iter] = funScatter2Grid2D(tempXY_refB(:,1),tempXY_refB(:,2),tempUV_B2A(:,1),sxy,0,xGrid,yGrid);
            [~,~,vGrid_B2A_refB_iter] = funScatter2Grid2D(tempXY_refB(:,1),tempXY_refB(:,2),tempUV_B2A(:,2),sxy,0,xGrid,yGrid);
        end
        
        [M,N] = size(uGrid_B2A_refB_iter);
        DMat = funDerivativeOp(M,N,mean(sxy));
        uVec = [uGrid_B2A_refB_iter(:),vGrid_B2A_refB_iter(:)]'; uVec=uVec(:);
        
        if iterNum==1
            vdualVec=0*uVec; mu=1; alphaList=[1e-2,1e-1,1e0,1e1,1e2,1e3];
            
            % close all; [coordinatesFEM_refB,elementsFEM_refB] = funMeshSetUp(xGrid,yGrid);
            for tempi = 1:length(alphaList)
                alpha = alphaList(tempi);
                uhatVec = (alpha*(DMat')*DMat + speye(2*M*N))\((uVec-vdualVec));
                % Plotdisp_show(uhatVec,coordinatesFEM_refB,elementsFEM_refB);
                Err1(tempi) = sqrt((uhatVec-uVec+vdualVec)'*(uhatVec-uVec+vdualVec));
                Err2(tempi) = sqrt((DMat*uhatVec)'*(DMat*uhatVec));
            end
            
            ErrSum = Err1/max(Err1) + Err2/max(Err2);
            [~,indexOfalpha] = min(ErrSum);
            try % Tune the best beta by a quadratic polynomial fitting
                [fitobj] = fit(log10(alphaList(indexOfalpha-1:1:indexOfalpha+1))',ErrSum(indexOfalpha-1:1:indexOfalpha+1),'poly2');
                p = coeffvalues(fitobj); alpha_best = 10^(-p(2)/2/p(1));
            catch, alpha_best = alphaList(indexOfalpha);
            end
        end
        
        
        uhatVec = (alpha_best*(DMat')*DMat + speye(2*M*N))\((uVec-vdualVec));
        
        vdualVec = vdualVec + uhatVec - uVec;
        
        uGrid_B2A_refB_iter = reshape( uhatVec(1:2:end), size(uGrid_B2A_refB_iter) );
        vGrid_B2A_refB_iter = reshape( uhatVec(2:2:end), size(uGrid_B2A_refB_iter) );
        
        % figure, plot(Err1/max(Err1),Err2/max(Err2),'o-');
        % figure, plot([1:5],ErrSum,'o-');
        
    end % END of "if gbSolver"
    
    
end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%% Store previous results
xyGrid_prevCurr = [xGrid(:), yGrid(:)];
uvGrid_B2A_refB_prevCurr = [uGrid_B2A_refB(:), vGrid_B2A_refB(:)];

%%%%% Find "parCoordB" in the raw deformed image %%%%%
Fx = scatteredInterpolant(xyGrid_prevCurr,uGrid_B2A_refB(:),'linear','linear');
tempu = Fx(parCoordBCurr);
Fy = scatteredInterpolant(xyGrid_prevCurr,vGrid_B2A_refB(:),'linear','linear');
tempv = Fy(parCoordBCurr);
parCoordBCurr = parCoordBCurr - [tempu(:),tempv(:)];

close all;
figure,plotCone2(xGrid(:), yGrid(:), uGrid_B2A_refB(:), vGrid_B2A_refB(:)); view(2);
set(gca,'fontsize',18); view(3); box on; axis equal; axis tight;  view(2);
title('Tracked disp field','fontweight','normal');


end



