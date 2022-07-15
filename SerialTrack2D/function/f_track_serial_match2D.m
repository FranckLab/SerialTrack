function [matches_A2B,u_B2A_curr_refB,track_A2B] = f_track_serial_match2D( parCoordA, parCoordB, varargin )
%FUNCTION matches = f_track_neightopo_match2(parCoordA,parCoordB,f_o_s,n_neighbors)
% Objective: Tracking based on topology similarity
% ---------------------------------------------------
%
%   INPUT:      parCoordA        - coordinates of particles in image A [n x 2]
%    
%               parCoordB        - coordinates of particles in image B [n x 2]
% 
%               f_o_s            - field of search [px]
%
%               n_neighbours     - number of neighbouring particles [integer]
%
%   OUTPUT:     matches_A2B      - list of indices of matching particles [m x 2]
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
f_o_s = 30;              % Size of search field [px]
n_neighborsMax = 25;     % Max # of neighboring particles
n_neighborsMin = 1;      % Min # of neighboring particles
locSolver = 1;           % Local solver: 1-topology-based feature; 2-histogram-based feature first and then topology-based feature;
gbSolver = 1;            % Global step solver: 1-moving least square fitting; 2-global regularization; 3-ADMM iterations
smoothness = 1e-2;       % Coefficient of regularization
outlrThres = 2;          % Threshold for removing outliers in TPT (Westerweel et al., Exp.Fluids, 2005)
maxIterNum = 20;         % Max ADMM iteration number
iterStopThres = 1e-3;    % ADMM iteration stopping threshold
usePrevResults = 0;      % Whether use previous results or not  
strain_f_o_s = 50;       % Size of virtual strain gauge [px]
strain_n_neighbors = 20; % # of neighboring particles used in strain gauge
MPTparaDistMissing = 2; % Distance threshold to check whether particle has a match or not [px]
gridxyROIRange.gridx = [min(parCoordA(:,1)),max(parCoordA(:,1))]; % ROI-gridx
gridxyROIRange.gridy = [min(parCoordA(:,2)),max(parCoordA(:,2))]; % ROI-gridy

% The user defined parameters
for tempi = 1:floor(length(varargin)/2)
     
    switch lower(varargin{2*tempi-1})
        case lower('f_o_s')
            f_o_s = varargin{2*tempi};
        case lower('n_neighborsMax')
            n_neighborsMax = varargin{2*tempi};
        case lower('n_neighborsMin')
            n_neighborsMin = varargin{2*tempi};
        case lower('locSolver')
            locSolver = varargin{2*tempi};    
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
        case lower('gridxyROIRange')
            gridxyROIRange = varargin{2*tempi};
        case lower('parCoordB_prev') % ### Only used in cumulative mode and using previous results
            parCoordB_prev = varargin{2*tempi};
        case lower('uv_B2A_prev') % ### Only used in cumulative mode and using previous results
            uv_B2A_prev = varargin{2*tempi};
        case lower('ImgSeqNum') % ### Only used in cumulative mode and using previous results
            ImgSeqNum = varargin{2*tempi};
        case lower('MPTparaDistMissing')
            MPTparaDistMissing = varargin{2*tempi};
        otherwise
    end
end
    
sxy = min([round(0.5*f_o_s),20])*[1,1]; % Grid size for regularization


    
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
    [tempu,tempv] = funInitGuess(parCoordB_prev,uv_B2A_prev,parCoordBCurr,ImgSeqNum);
    u_B2A_curr_refB = u_B2A_curr_refB + [tempu,tempv];
    parCoordBCurr = parCoordBCurr + [tempu,tempv];
end


%% %%%%%%%%%%%%%%%%%%%%%% ADMM iterations %%%%%%%%%%%%%%%%%%%%% 
while iterNum < maxIterNum
    
    % close all;
    % # of neighboring particles: >1 (topology) => 1 (nearest neighboring search)
    n_neighbors = round(n_neighborsMin + exp(-0.5*iterNum)*(n_neighborsMax-n_neighborsMin)); 
    % current iteration number
    iterNum = iterNum+1; disp(['------ Iter #',num2str(iterNum),' ------']); 
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Local step: neighbor topology match %%%%%
    matches_A2B = []; 
    while isempty(matches_A2B) && n_neighborsMax < length(parNotMissingIndA)
        
        if locSolver == 1 && n_neighbors > 2 % topology-based feature
            matches_A2B = f_track_neightopo_match2( parCoordA(parNotMissingIndA,:), parCoordBCurr(parNotMissingIndBCurr,:), f_o_s, n_neighbors );
        elseif locSolver == 2 && n_neighbors > 2 % histogram-based feature
            gauss_interp = 1; matches_A2B = f_track_hist_match( parCoordA(parNotMissingIndA,:), parCoordBCurr(parNotMissingIndBCurr,:), f_o_s, n_neighbors, gauss_interp );
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
    if isempty(matches_A2B)==1, disp('No matches found!'); break; end
    
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Global step: Compute kinematically compatible displacement %%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Method I: local moving least squares %%%%%
    if gbSolver==1
        
        [XY_B2A_refB,U_B2A_refB,F_B2A_refB] = funCompDefGrad2(-u_A2B, parCoordBCurr(track_A2B(track_A2B>0),:), f_o_s, n_neighbors);
        
        [row,~] = find(isnan(U_B2A_refB(1:2:end)) == 1);        % find nans
        XY_B2A_refB(row,:) = [];                                % remove nans
        U_B2A_refB([2*row-1; 2*row]) = [];                      % remove nans
        F_B2A_refB([4*row-3; 4*row-2; 4*row-1; 4*row]) = [];    % remove nans
        
        Fx = scatteredInterpolant(XY_B2A_refB,U_B2A_refB(1:2:end),'linear','linear');
        tempu = Fx(parCoordBCurr); % Interpolate temporary x-displacement
        Fy = scatteredInterpolant(XY_B2A_refB,U_B2A_refB(2:2:end),'linear','linear');
        tempv = Fy(parCoordBCurr); % Interpolate temporary y-displacement
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Method II: global regularization %%%%%    
    elseif gbSolver==2
        
        tempUV_B2A = -u_A2B;
        tempXY_refB = parCoordBCurr(track_A2B(track_A2B>0),:);
        try
            [xGrid_refB,yGrid_refB,uGrid_B2A_refB_iter]=funScatter2Grid2D(tempXY_refB(:,1),tempXY_refB(:,2),tempUV_B2A(:,1),sxy,smoothness);
            [~,~,vGrid_B2A_refB_iter]=funScatter2Grid2D(tempXY_refB(:,1),tempXY_refB(:,2),tempUV_B2A(:,2),sxy,smoothness);
        catch
            [xGrid_refB,yGrid_refB,uGrid_B2A_refB_iter]=funScatter2Grid2D(tempXY_refB(:,1),tempXY_refB(:,2),tempUV_B2A(:,1),sxy,0);
            [~,~,vGrid_B2A_refB_iter]=funScatter2Grid2D(tempXY_refB(:,1),tempXY_refB(:,2),tempUV_B2A(:,2),sxy,0);
        end
        Fx = scatteredInterpolant([xGrid_refB(:),yGrid_refB(:)],uGrid_B2A_refB_iter(:),'linear','linear');
        tempu = Fx(parCoordBCurr); % Interpolate temporary x-displacement
        Fy = scatteredInterpolant([xGrid_refB(:),yGrid_refB(:)],vGrid_B2A_refB_iter(:),'linear','linear');
        tempv = Fy(parCoordBCurr); % Interpolate temporary y-displacement
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Method III: augmented Lagrangian regularization %%%%%
    elseif gbSolver==3
        
        tempUV_B2A = -u_A2B;
        tempXY_refB = parCoordBCurr(track_A2B(track_A2B>0),:);
        if iterNum==1
            [xGrid_refB,yGrid_refB] = ndgrid(gridxyROIRange.gridx(1) : sxy(1) : gridxyROIRange.gridx(2), ...
                                             gridxyROIRange.gridy(1) : sxy(2) : gridxyROIRange.gridy(2));
        end
        try
            [~,~,uGrid_B2A_refB_iter]=funScatter2Grid2D(tempXY_refB(:,1),tempXY_refB(:,2),tempUV_B2A(:,1),sxy,smoothness,xGrid_refB,yGrid_refB);
            [~,~,vGrid_B2A_refB_iter]=funScatter2Grid2D(tempXY_refB(:,1),tempXY_refB(:,2),tempUV_B2A(:,2),sxy,smoothness,xGrid_refB,yGrid_refB);
        catch
            [~,~,uGrid_B2A_refB_iter]=funScatter2Grid2D(tempXY_refB(:,1),tempXY_refB(:,2),tempUV_B2A(:,1),sxy,0,xGrid_refB,yGrid_refB);
            [~,~,vGrid_B2A_refB_iter]=funScatter2Grid2D(tempXY_refB(:,1),tempXY_refB(:,2),tempUV_B2A(:,2),sxy,0,xGrid_refB,yGrid_refB);
        end
    
        [M,N] = size(uGrid_B2A_refB_iter); % Grid size [M,N]
        DMat = funDerivativeOp(M,N,mean(sxy)); % Build a finite different operator such that {F}={D}{U}
        uVec = [uGrid_B2A_refB_iter(:),vGrid_B2A_refB_iter(:)]'; uVec=uVec(:);
        
        %%%%% Tune parameter alpha by L-curve method %%%%%
        if iterNum==1 
            vdualVec=0*uVec; mu=1; alphaList=[1e-2,1e-1,1e0,1e1,1e2,1e3];
            for tempi = 1:length(alphaList)
                alpha = alphaList(tempi);
                uhatVec = (alpha*(DMat')*DMat + speye(2*M*N))\((uVec-vdualVec));
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
        uhatVec = (alpha_best*(DMat')*DMat + speye(2*M*N))\((uVec-vdualVec));
        
        %%%%% Update dual variable
        vdualVec = vdualVec + uhatVec - uVec;
        
        %%%%% Interpolate to scatterred data points %%%%% 
        uGrid_B2A_refB_iter = reshape( uhatVec(1:2:end), size(uGrid_B2A_refB_iter) );
        vGrid_B2A_refB_iter = reshape( uhatVec(2:2:end), size(uGrid_B2A_refB_iter) );
        
        Fx = scatteredInterpolant([xGrid_refB(:),yGrid_refB(:)],uGrid_B2A_refB_iter(:),'linear','linear');
        tempu = Fx(parCoordBCurr); % Interpolate temporary x-displacement
        Fy = scatteredInterpolant([xGrid_refB(:),yGrid_refB(:)],vGrid_B2A_refB_iter(:),'linear','linear');
        tempv = Fy(parCoordBCurr); % Interpolate temporary y-displacement
        
            
    end % END of if gbSolver==0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Check convergence %%%%%
    %%%%% Norm of disp update %%%%%
    u_B2ACurr_updateNorm = sqrt( ( norm(tempu)^2+norm(tempv)^2 )/length(tempu) );
    disp(['Disp update norm: ',num2str(u_B2ACurr_updateNorm)]);
    
    %%%%% Do at most five iterations that "matchRatio==1" %%%%%
    if matchRatio > 0.999
        matchRatioEqualsOneTime = matchRatioEqualsOneTime+1;
    end
    
    %%%%% Stopping criterion %%%%%
    if u_B2ACurr_updateNorm < sqrt(2)*iterStopThres || matchRatioEqualsOneTime>5
        disp(['----- Converged! ------']); break
    
    else  
        %%%%% Warp back parCoordBCurr %%%%%
        u_B2A_curr_refB = u_B2A_curr_refB + [tempu,tempv ];
        parCoordBCurr = parCoordBCurr + [tempu,tempv ];
        
        %%%%% Update f_o_s %%%%%
        tempu_Quantile = quantile(tempu,[0.25,0.5,0.75]);
        tempv_Quantile = quantile(tempv,[0.25,0.5,0.75]);
        f_o_s = max( [ 30; tempu_Quantile(2)+0.5*(tempu_Quantile(3)-tempu_Quantile(1));
                       tempv_Quantile(2)+0.5*(tempv_Quantile(3)-tempv_Quantile(1))]);
           
        %%%%% Remove non-paired particles because of particle missing or merging %%%%%           
        if n_neighbors < 4 % n_neighborsMax/2            
            neighborInd_BCurrA = knnsearch(parCoordBCurr(:,1:2),parCoordA(:,1:2),'K',1); % Find pts in BCurr near pts in ACurr
            dist_BCurrA = sqrt( sum((parCoordBCurr(neighborInd_BCurrA,1:2) -  parCoordA).^2,2) ); % figure, h = histogram(dist_BCurrA);
            [parNotMissingIndA,~] = find(dist_BCurrA < max([2, MPTparaDistMissing])); % Find particles not missing in Particle A
            
            neighborInd_ABCurr = knnsearch(parCoordA(:,1:2),parCoordBCurr(:,1:2),'K',1); % Find pts in ACurr near pts in BCurr
            dist_ABCurr = sqrt( sum((parCoordA(neighborInd_ABCurr,1:2) -  parCoordBCurr).^2,2) );
            [parNotMissingIndBCurr,~] = find(dist_ABCurr < max([2, MPTparaDistMissing]));  % Find particles not missing in Particle BCurr
        end           
        
        %%%%% Plot detected particles %%%%%
        % close all; figure,
        % hold on; plot( parCoordA(:,1), parCoordA(:,2), 'bo');
        % hold on; plot( parCoordA(parNotMissingIndA,1), parCoordA(parNotMissingIndA,2), 'b.');
        % view(2); box on; axis equal; axis tight; set(gca,'fontsize',18);
        % title('Detected particles in ref image','fontweight','normal');
        % 
        % % figure,
        % hold on; plot(parCoordBCurr(:,1), parCoordBCurr(:,2),'ro');
        % hold on; plot(parCoordBCurr(parNotMissingIndBCurr,1), parCoordBCurr(parNotMissingIndBCurr,2),'r.');
        % view(2); box on; axis equal; axis tight; set(gca,'fontsize',18);
        % title('Detected particles in defor image','fontweight','normal');
        % pause;

    end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%% SerialTrack: Done! %%%%%%'); 
timeCost = toc; toc
fprintf('\n'); 












