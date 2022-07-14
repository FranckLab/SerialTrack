function [matches_A2B,u_B2A_curr_refB,track_A2B] = f_track_serial_match3D( parCoordA, parCoordB, varargin )
%FUNCTION matches = f_track_serial_match3D(parCoordA,parCoordB,f_o_s,n_neighbors)
% Objective: Tracking based on topology similarity
% ---------------------------------------------------
%
%   INPUT:      parCoordA        - coordinates of particles in image A [n x 3]
%
%               parCoordB        - coordinates of particles in image B [n x 3]
%
%               f_o_s         - field of search [px]
%
%               n_neighbours  - number of neighbouring particles [integer]
%
%   OUTPUT:     matches_A2B       - list of indices of matching particles [m x 3]
%
% ---------------------------------------------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%
% ---------------------------------------------------
% Author: Jin Yang
% Contact and support: jyang526@wisc.edu
% Date: 2020.12.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% %%%%%%%%% Parameter inputs %%%%%%%%%%%
% Default parameters
f_o_s = 60;              % Size of search field [px]
n_neighborsMax = 25;     % Max # of neighboring particles
n_neighborsMin = 1;      % Min # of neighboring particles
gbSolver = 1;            % Global step solver: 1-moving least square fitting; 2-global regularization; 3-ADMM iterations
smoothness = 1e-1;       % Coefficient of regularization
outlrThres = 5;          % Threshold for removing outliers in TPT (Westerweel et al., Exp.Fluids, 2005)
maxIterNum = 20;         % Max ADMM iteration number
iterStopThres = 1e-2;    % ADMM iteration stopping threshold
usePrevResults = 0;      % Whether use previous results or not
strain_f_o_s = 60;       % Size of virtual strain gauge [px]
strain_n_neighbors = 20; % # of neighboring particles used in strain gauge
MPTParaDistMissing = 2;  % Distance threshold to check whether particle has a match or not [px]
gridxyzROIRange.gridx = [min(parCoordA(:,1)),max(parCoordA(:,1))]; % ROI-gridx
gridxyzROIRange.gridy = [min(parCoordA(:,2)),max(parCoordA(:,2))]; % ROI-gridy
gridxyzROIRange.gridz = [min(parCoordA(:,3)),max(parCoordA(:,3))]; % ROI-gridz

% The user defined parameters
for tempi = 1:floor(length(varargin)/2)
    
    switch lower(varargin{2*tempi-1})
        case lower('f_o_s')
            f_o_s = varargin{2*tempi};
        case lower('n_neighborsMax')
            n_neighborsMax = varargin{2*tempi};
        case lower('n_neighborsMin')
            n_neighborsMin = varargin{2*tempi};
        case lower('gbSolver')
            gbSolver = varargin{2*tempi};
        case lower('smoothness')
            smoothness = varargin{2*tempi};
        case lower('outlrThres')
            outlrThres = varargin{2*tempi};
        case lower('maxIterNum')
            maxIterNum = varargin{2*tempi};
        case lower('iterStopThres')
            iterStopThres = varargin{2*tempi};
        case lower('usePrevResults')
            usePrevResults = varargin{2*tempi};
        case lower('strain_f_o_s')
            strain_f_o_s = varargin{2*tempi};
        case lower('strain_n_neighbors')
            strain_n_neighbors = varargin{2*tempi};
        case lower('gridxyzROIRange')
            gridxyzROIRange = varargin{2*tempi};
        case lower('parCoordB_prev') % ### Only used in cumulative mode and using previous results
            parCoordB_prev = varargin{2*tempi};
        case lower('uvw_B2A_prev') % ### Only used in cumulative mode and using previous results
            uvw_B2A_prev = varargin{2*tempi};
        case lower('ImgSeqNum') % ### Only used in cumulative mode and using previous results
            ImgSeqNum = varargin{2*tempi};
        case lower('MPTParaDistMissing')
            MPTParaDistMissing = varargin{2*tempi};
        otherwise
    end
end

sxyz = min([round(0.5*f_o_s),20])*[1,1,1]; % Grid size for regularization



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Initialize iteration %%%%%
tic; iterNum = 0; % Initialize ADMM iteration number
parCoordBCurr = parCoordB; % Current particle B coordinates
u_B2A_curr_refB = 0*parCoordBCurr; % Current displacements from image B to image A
parNotMissingIndA = [1:size(parCoordA,1)]'; % Indices of particles in image A that have matches in image B
parNotMissingIndBCurr = [1:size(parCoordBCurr,1)]'; % Indices of particles in image B that have matches in image A
matchRatioEqualsOneTime = 0; % How many times the match raio equals one


%%%%% Use previous results for a new frame %%%%%
if ImgSeqNum>2 && usePrevResults
    % disp('Use previous results as an initial estimate');
    [tempu,tempv,tempw] = funInitGuess3(parCoordB_prev,uvw_B2A_prev,parCoordBCurr,ImgSeqNum);
    u_B2A_curr_refB = u_B2A_curr_refB + [tempu,tempv,tempw];
    parCoordBCurr = parCoordBCurr + [tempu,tempv,tempw];
end


%% %%%%%%%%%%%%%%%%%%%%%% ADMM iterations %%%%%%%%%%%%%%%%%%%%%
while iterNum < maxIterNum
    
    close all;
    % # of neighboring particles: >1 (topology) => 1 (nearest neighboring search)
    n_neighbors = round(n_neighborsMin + exp(-0.5*iterNum)*(n_neighborsMax-n_neighborsMin));
    % current iteration number
    iterNum = iterNum+1; disp(['------ Iter #',num2str(iterNum),' ------']);
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Local step: neighbor topology match %%%%%
    matches_A2B = [];
    while isempty(matches_A2B) && n_neighborsMax < length(parNotMissingIndA)
        
        if n_neighbors > 2
            matches_A2B = f_track_neightopo_match3( parCoordA(parNotMissingIndA,:), parCoordBCurr(parNotMissingIndBCurr,:), f_o_s, n_neighbors );
            % matches_A2B = f_track_hist_match( parCoordA(parNotMissingIndA,:), parCoordBCurr(parNotMissingIndBCurr,:), f_o_s, n_neighbors, gauss_interp );
        else
            matches_A2B = f_track_nearest_neighbour3( parCoordA(parNotMissingIndA,:), parCoordBCurr(parNotMissingIndBCurr,:), f_o_s );
        end
        if isempty(matches_A2B) == 1
            n_neighborsMax = round(n_neighborsMax + 5);
        else
            matches_A2B = [parNotMissingIndA(matches_A2B(:,1)), parNotMissingIndBCurr(matches_A2B(:,2))];
            [track_A2B, u_A2B] = funCompDisp3(parCoordA, parCoordBCurr, matches_A2B, outlrThres);
            matchRatio = size(matches_A2B,1) / length(parNotMissingIndA);
            disp( ['Tracking ratio: ', num2str( size(matches_A2B,1) ),'/', num2str(length(parNotMissingIndA)), ' = ', num2str(matchRatio) ]);
        end
    end
    if isempty(matches_A2B)==1, disp('Wrong parameters!'); break; end
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Global step: Compute kinematically compatible displacement %%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Method I: local moving least squares %%%%%
    if gbSolver==1
        
        [XYZ_B2A_refB,U_B2A_refB,F_B2A_refB] = funCompDefGrad3(-u_A2B, parCoordBCurr(track_A2B(track_A2B>0),:), f_o_s, n_neighbors);
        
        [row,~] = find(isnan(U_B2A_refB(1:3:end)) == 1);         % find nans
        XYZ_B2A_refB(row,:) = [];                                % remove nans
        U_B2A_refB([3*row-2; 3*row-1; 3*row]) = [];              % remove nans
        F_B2A_refB([9*row-8; 9*row-7; 9*row-6; 9*row-5;  ...
            9*row-4; 9*row-3; 9*row-2; 9*row-1; 9*row]) = [];    % remove nans
        
        Fx = scatteredInterpolant(XYZ_B2A_refB,U_B2A_refB(1:3:end),'linear','linear');
        tempu = Fx(parCoordBCurr); % Interpolate temporary x-displacement
        Fy = scatteredInterpolant(XYZ_B2A_refB,U_B2A_refB(2:3:end),'linear','linear');
        tempv = Fy(parCoordBCurr); % Interpolate temporary y-displacement
        Fz = scatteredInterpolant(XYZ_B2A_refB,U_B2A_refB(3:3:end),'linear','linear');
        tempw = Fz(parCoordBCurr); % Interpolate temporary z-displacement
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Method II: global regularization %%%%%
    elseif gbSolver==2
        
        tempUVW_B2A = -u_A2B;
        tempXYZ_refB = parCoordBCurr(track_A2B(track_A2B>0),:);
        try
            try
                [xGrid_refB,yGrid_refB,zGrid_refB,uGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,1),sxyz,smoothness);
                [~,~,~,vGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,2),sxyz,smoothness);
                [~,~,~,wGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,3),sxyz,smoothness);
            catch
                [xGrid_refB,yGrid_refB,zGrid_refB,uGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,1),sxyz,0);
                [~,~,~,vGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,2),sxyz,0);
                [~,~,~,wGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,3),sxyz,0);
            end
        catch
            [xGrid_refB,yGrid_refB,zGrid_refB,uGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,1),[1,1,1],0);
            [~,~,~,vGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,2),[1,1,1],0);
            [~,~,~,wGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,3),[1,1,1],0);
        end
        Fx = scatteredInterpolant([xGrid_refB(:),yGrid_refB(:),zGrid_refB(:)],uGrid_B2A_refB_iter(:),'linear','linear');
        tempu = Fx(parCoordBCurr);
        Fy = scatteredInterpolant([xGrid_refB(:),yGrid_refB(:),zGrid_refB(:)],vGrid_B2A_refB_iter(:),'linear','linear');
        tempv = Fy(parCoordBCurr);
        Fz = scatteredInterpolant([xGrid_refB(:),yGrid_refB(:),zGrid_refB(:)],wGrid_B2A_refB_iter(:),'linear','linear');
        tempw = Fz(parCoordBCurr);
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Method III: augmented Lagrangian regularization %%%%%
    elseif gbSolver==3
        
        tempUVW_B2A = -u_A2B;
        tempXYZ_refB = parCoordBCurr(track_A2B(track_A2B>0),:);
        if iterNum==1
            [xGrid_refB,yGrid_refB,zGrid_refB] = ndgrid(TPTpara.gridxyzROIRange.gridx(1) : sxyz(1) : TPTpara.gridxyzROIRange.gridx(2), ...
                TPTpara.gridxyzROIRange.gridy(1) : sxyz(2) : TPTpara.gridxyzROIRange.gridy(2), ...
                TPTpara.gridxyzROIRange.gridz(1) : sxyz(3) : TPTpara.gridxyzROIRange.gridz(2) );
        end
        try
            try
                [~,~,~,uGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,1),sxyz,smoothness,xGrid_refB,yGrid_refB,zGrid_refB);
                [~,~,~,vGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,2),sxyz,smoothness,xGrid_refB,yGrid_refB,zGrid_refB);
                [~,~,~,wGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,3),sxyz,smoothness,xGrid_refB,yGrid_refB,zGrid_refB);
            catch
                [~,~,~,uGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,1),sxyz,0,xGrid_refB,yGrid_refB,zGrid_refB);
                [~,~,~,vGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,2),sxyz,0,xGrid_refB,yGrid_refB,zGrid_refB);
                [~,~,~,wGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,3),sxyz,0,xGrid_refB,yGrid_refB,zGrid_refB);
            end
        catch
            [~,~,~,uGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,1),[1,1,1],0,xGrid_refB,yGrid_refB,zGrid_refB);
            [~,~,~,vGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,2),[1,1,1],0,xGrid_refB,yGrid_refB,zGrid_refB);
            [~,~,~,wGrid_B2A_refB_iter]=funScatter2Grid3D(tempXYZ_refB(:,1),tempXYZ_refB(:,2),tempXYZ_refB(:,3),tempUVW_B2A(:,3),[1,1,1],0,xGrid_refB,yGrid_refB,zGrid_refB);
        end
        
        [M,N,L] = size(uGrid_B2A_refB_iter); % Grid size [M,N,L]
        DMat = funDerivativeOp3(M,N,L,mean(sxyz)*[1,1,1]); % Build a finite different operator such that {F}={D}{U}
        uVec = [uGrid_B2A_refB_iter(:),vGrid_B2A_refB_iter(:),wGrid_B2A_refB_iter(:)]'; uVec=uVec(:);
        
        %%%%% Tune parameter alpha by L-curve method %%%%%
        if iterNum==1
            vdualVec=0*uVec; mu=1; alphaList=[1e-2,1e-1,1e0,1e1,1e2,1e3];
            for tempi = 1:length(alphaList)
                alpha = alphaList(tempi);
                uhatVec = (alpha*(DMat')*DMat + speye(3*M*N*L))\((uVec-vdualVec));
                Err1(tempi) = sqrt((uhatVec-uVec+vdualVec)'*(uhatVec-uVec+vdualVec));
                Err2(tempi) = sqrt((DMat*uhatVec)'*(DMat*uhatVec));
            end
            ErrSum = Err1/max(Err1) + Err2/max(Err2);
            [~,indexOfalpha] = min(ErrSum);
            try % Find the best value of alpha by a quadratic polynomial fitting
                [fitobj] = fit(log10(alphaList(indexOfalpha-1:1:indexOfalpha+1))',ErrSum(indexOfalpha-1:1:indexOfalpha+1),'poly2');
                p = coeffvalues(fitobj); alpha_best = 10^(-p(2)/2/p(1));
            catch, alpha_best = alphaList(indexOfalpha);
            end
        end
        
        %%%%% Resolve global step with tuned best alpha %%%%%
        uhatVec = (alpha_best*(DMat')*DMat + speye(3*M*N*L))\((uVec-vdualVec));
        
        %%%%% Update dual variable
        vdualVec = vdualVec + uhatVec - uVec;
        
        %%%%% Interpolate to scatterred data points %%%%%
        uGrid_B2A_refB_iter = reshape( uhatVec(1:3:end), size(uGrid_B2A_refB_iter) );
        vGrid_B2A_refB_iter = reshape( uhatVec(2:3:end), size(uGrid_B2A_refB_iter) );
        wGrid_B2A_refB_iter = reshape( uhatVec(3:3:end), size(uGrid_B2A_refB_iter) );
        
        Fx = scatteredInterpolant([xGrid_refB(:),yGrid_refB(:),zGrid_refB(:)],uGrid_B2A_refB_iter(:),'linear','linear');
        tempu = Fx(parCoordBCurr);
        Fy = scatteredInterpolant([xGrid_refB(:),yGrid_refB(:),zGrid_refB(:)],vGrid_B2A_refB_iter(:),'linear','linear');
        tempv = Fy(parCoordBCurr);
        Fz = scatteredInterpolant([xGrid_refB(:),yGrid_refB(:),zGrid_refB(:)],wGrid_B2A_refB_iter(:),'linear','linear');
        tempw = Fz(parCoordBCurr);
        
        
    end % END of if gbSolver==0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Check convergence %%%%%
    %%%%% Norm of disp update %%%%%
    u_B2ACurr_updateNorm = sqrt( ( norm(tempu)^2+norm(tempv)^2+norm(tempw)^2 )/length(tempu) );
    disp(['Disp update norm: ',num2str(u_B2ACurr_updateNorm)]);
    
    %%%%% Do at most five iterations that "matchRatio==1" %%%%%
    if matchRatio > 0.999
        matchRatioEqualsOneTime = matchRatioEqualsOneTime+1;
    end
    
    %%%%% Stopping criterion %%%%%
    if u_B2ACurr_updateNorm < sqrt(3)*iterStopThres || matchRatioEqualsOneTime>5
        disp(['----- Converged! ------']); break
        
    else
        %%%%% Warp back parCoordBCurr %%%%%
        u_B2A_curr_refB = u_B2A_curr_refB + [tempu,tempv,tempw];
        parCoordBCurr = parCoordBCurr + [tempu,tempv,tempw];
        
        %%%%% Update f_o_s %%%%%
        tempu_Quantile = quantile(tempu,[0.25,0.5,0.75]);
        tempv_Quantile = quantile(tempv,[0.25,0.5,0.75]);
        tempw_Quantile = quantile(tempw,[0.25,0.5,0.75]);
        f_o_s = max( [ 60; tempu_Quantile(2)+0.5*(tempu_Quantile(3)-tempu_Quantile(1));
                           tempv_Quantile(2)+0.5*(tempv_Quantile(3)-tempv_Quantile(1));
                           tempw_Quantile(2)+0.5*(tempw_Quantile(3)-tempw_Quantile(1))]);
        
        %%%%% Remove non-paired particles because of particle missing or merging %%%%%
        if n_neighbors < 4 %  n_neighborsMax/2
            neighborInd_BCurrA = knnsearch(parCoordBCurr(:,1:3),parCoordA(:,1:3),'K',1); % Find pts in BCurr near pts in ACurr
            dist_BCurrA = sqrt( sum((parCoordBCurr(neighborInd_BCurrA,1:3) -  parCoordA).^2,2) ); % figure, h = histogram(dist_BCurrA);
            [parNotMissingIndA,~] = find(dist_BCurrA < max([2, MPTParaDistMissing])); % Find particles not missing in Particle A
            
            neighborInd_ABCurr = knnsearch(parCoordA(:,1:3),parCoordBCurr(:,1:3),'K',1); % Find pts in ACurr near pts in BCurr
            dist_ABCurr = sqrt( sum((parCoordA(neighborInd_ABCurr,1:3) - parCoordBCurr).^2,2) );
            
            [parNotMissingIndBCurr,~] = find(dist_ABCurr < max([2, MPTParaDistMissing]));  % Find particles not missing in Particle BCurr
        end
        
        
    end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%% SerialTrack: Done! %%%%%%');
timeCost = toc; toc
fprintf('\n');












