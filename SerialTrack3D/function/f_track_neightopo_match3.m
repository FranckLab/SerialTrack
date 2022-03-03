function [matches] = f_track_neightopo_match3( part_A, part_B, f_o_s, n_neighbors)
%FUNCTION matches = f_track_neightopo_match3(part_A,part_B,f_o_s,n_neighbors)
% Objective: Tracking based on topology similarity
% ---------------------------------------------------
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
% ---------------------------------------------------
%
%   INPUT:      part_A        - coordinates of particles in image A [n x 3]
%
%               part_B        - coordinates of particles in image B [n x 3]
%
%               f_o_s         - field of search [px]
%
%               n_neighbours  - number of neighbouring particles [integer]
%
%   OUTPUT:     matches       - list of indices of matching particles [m x 2]
%
% ---------------------------------------------------
% References
% [1] T Janke, R Schwarze, K Bauer. Part2Track: A MATLAB package for double
%     frame and time resolved Particle Tracking Velocimetry. 11, 100413, SoftwareX (2020).
% ---------------------------------------------------
% Author: Jin Yang
% Contact and support: jyang526@wisc.edu
% Date: 2020.12.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%% Example values %%%%%
% part_A = partA;
% part_B = partB;
% n_neighbors = 20;
% f_o_s = 30;
% f_o_s = Inf;

%%
% Calculate neighbourhood indices for each particle in image A and image B
for tempAorB = 1:2
    
    if tempAorB==1 % Image A
        neighborInd = knnsearch(part_A(:,1:3),part_A(:,1:3),'K',n_neighbors+1);
        part_curr = part_A;
    elseif tempAorB==2 % Image B
        neighborInd = knnsearch(part_B(:,1:3),part_B(:,1:3),'K',n_neighbors+1);
        part_curr = part_B;
    end
    
    featureVec_r = zeros(size(part_curr,1),n_neighbors);
    featureVec_phi = zeros(size(part_curr,1),n_neighbors);
    featureVec_theta = zeros(size(part_curr,1),n_neighbors);
    
    for parInd = 1:size(part_curr,1) % For each particle in part_curr
        
        % Grab neighboring particles
        parNeighborIndCoord = part_curr(neighborInd(parInd,:),:);
        
        % Compute distance btw neighboring particles to particle "parInd"
        dist_xyz = parNeighborIndCoord(2:end,:) - repmat(part_curr(parInd,:), n_neighbors, 1);
        
        % Define rotation invariant (RI) coordinate system
        ex = dist_xyz(1,:) /  sqrt(sum(dist_xyz(1,:).^2)) ; % ex is along the nearest neighbor
        ez = [det([dist_xyz(1,2:3); dist_xyz(2,2:3)]), -det([dist_xyz(1,[1,3]); dist_xyz(2,[1,3])]), det([dist_xyz(1,1:2); dist_xyz(2,1:2)])];
        ez = ez / sqrt(sum(ez.^2)); % ez is perpenditular to r1 and r2;
        if (ez*dist_xyz(3,:)'<0), ez=-ez; end % third nearest neighbor is in the ez positive direction
        ey = [det([ez(2:3);ex(2:3)]), -det([ez([1,3]);ex([1,3])]), det([ez(1:2);ex(1:2)])];
        ey = ey / sqrt(sum(ey.^2));
        
        % Rotation invariant (RI) transformation
        dist_xyz = dist_xyz*[ex', ey', ez'];
        
        % Compute topology feature {r, phi, theta}
        r = sqrt( sum(dist_xyz.^2,2) ) ; % abs distance 'r' is already in order
        phi =  atan2( dist_xyz(:,2) , dist_xyz(:,1) ) ; % phi angle in xy-plane
        theta = atan2( dist_xyz(:,3), sqrt( sum(dist_xyz(:,1:2).^2,2) )  ); % polar angle in z-axis
        
        [ dist_phi_sort, dist_phi_sort_ind ] = sort(phi); % sort angle phi
        tempInd = find( dist_phi_sort_ind == 1 ); % find the index of r_min
        tempInd2 = [ tempInd : length(dist_phi_sort),  1:(tempInd-1) ]; % reorder index starting from the r_min
        
        r2 = r( dist_phi_sort_ind( tempInd2 ) ); % Reorder distance 'r' starting from r_min
        phi2 = phi( dist_phi_sort_ind( tempInd2 ) ); % Reorder angle 'phi' starting from r_min
        phiDiff2 = phi2([2:end,1]) - phi2;  % Compute angle between two r's
        phiDiff2(phiDiff2 < 0) = phiDiff2(phiDiff2 < 0) + 2*pi; % Compensate 2*pi to make sure it's always positive
        theta2 = theta( dist_phi_sort_ind( tempInd2 )  ); % Reordered angle in z-axis
        
        featureVec_r(parInd,:) = [r2(:)']; % /min(r2); % Normalize r_min if needed
        featureVec_phi(parInd,:) = [phiDiff2(:)']; % Make sure the angle difference is a row array
        featureVec_theta(parInd,:) = [theta2(:)']; % Make sure the angle difference is a row array
        
    end
    
    % Store feature vectors
    if tempAorB==1 % Image A
        featureVec_r_A = featureVec_r; % feature vector in "r" dimension
        featureVec_phi_A = featureVec_phi; % feature vector in "phi" dim.,  sum of each row of featureVec_phi_A is 2*pi
        featureVec_theta_A = featureVec_theta;
    elseif tempAorB==2  % Image B
        featureVec_r_B = featureVec_r; % feature vector in "r" dimension
        featureVec_phi_B = featureVec_phi; % feature vector in "phi" dim.,
        featureVec_theta_B = featureVec_theta;
    end
    
    
end


%% Find the most similar particle in part_B for each particle in part_A
waitbarOrNot = 0; % Waitbar to be used to check the code progres
if waitbarOrNot==1
    hbar = parfor_progressbar(size(part_A,1),'Computing...');  % create the progress bar in parallel computing "parfor"
end
clear matches;
  
parfor parInd = 1:size(part_A,1) % For each particle in part_A
    
    % if waitbarOrNot==1
    %    hbar.iterate(1);
    % end
    
    %%%%% Only use neighbors within f_o_s %%%%%
    if f_o_s < Inf
        neighborInd_InBNearA = knnsearch(part_B(:,1:3),part_A(parInd,1:3),'K', n_neighbors);
        dist_neigh_InBNearA = sqrt(sum( (repmat(part_A(parInd,:),length(neighborInd_InBNearA),1) - part_B(neighborInd_InBNearA,:)).^2, 2));
        neighborInd_InBNearA = neighborInd_InBNearA(dist_neigh_InBNearA < sqrt(3)*f_o_s);
    else
        neighborInd_InBNearA = 1:size(part_B,1); % Search the whole field
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(neighborInd_InBNearA) == 0
        
        feaVec_r_A1 = featureVec_r_A(parInd,:);
        feaVec_phi_A1 = featureVec_phi_A(parInd,:);
        feaVec_theta_A1 = featureVec_theta_A(parInd,:);
        
        feaVec_r_B = featureVec_r_B(neighborInd_InBNearA(1:length(neighborInd_InBNearA)),:);
        feaVec_phi_B = featureVec_phi_B(neighborInd_InBNearA(1:length(neighborInd_InBNearA)),:);
        feaVec_theta_B = featureVec_theta_B(neighborInd_InBNearA(1:length(neighborInd_InBNearA)),:);
        
        % LSQ difference between each topology feature
        corr_r  = sum( (repmat(feaVec_r_A1,length(neighborInd_InBNearA),1)-feaVec_r_B).^2 , 2);
        corr_phi = sum( (repmat(feaVec_phi_A1,length(neighborInd_InBNearA),1)-feaVec_phi_B).^2 , 2);
        corr_theta = sum( (repmat(feaVec_theta_A1,length(neighborInd_InBNearA),1)-feaVec_theta_B).^2 , 2);
        
        [~,ind_min_corr_r] = min(corr_r);
        [~,ind_min_corr_phi] = min(corr_phi);
        [~,ind_min_corr_theta] = min(corr_theta);
        
        %%%%%% All three features (r,phi,theta) reach minimum %%%%%
        if (ind_min_corr_r==ind_min_corr_phi) && (ind_min_corr_phi==ind_min_corr_theta)
             
            % discard displacements should be smaller than the field of search
            if sqrt( sum( ( part_A(parInd,:) - part_B(neighborInd_InBNearA(ind_min_corr_r),:) ).^2 ) ) < f_o_s
                matches{parInd} = [parInd, neighborInd_InBNearA(ind_min_corr_r)]; % save matches
            end
            
            %%%%% Here I comment sometimes there can be a relaxed criterion %%%%% 
            %%%%%% Two of the three features reach minimum %%%%%
            %     if (ind_min_corr_r==ind_min_corr_phi)
            %
            %         % if .. < f_o_s % discard displacements larger field of search
            %         matches{parInd} = [parInd, neighborInd_AB(ind_min_corr_r)]; % save matches
            %
            %     elseif (ind_min_corr_phi==ind_min_corr_theta)
            %         matches{parInd} = [parInd, neighborInd_AB(ind_min_corr_phi)]; % save matches
            %
            %     elseif (ind_min_corr_r==ind_min_corr_theta)
            %         matches{parInd} = [parInd, neighborInd_AB(ind_min_corr_r)]; % save matches
            %
            %     end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
        
    end % END of if isempty(neighborInd_BA) == 0
    
    
end


% Check if there are found matches
if exist('matches','var')
    matches = cell2mat(matches');
else
    matches = [];
end

if waitbarOrNot==1, close(hbar); end % Close hbar if there was a waitbar.







