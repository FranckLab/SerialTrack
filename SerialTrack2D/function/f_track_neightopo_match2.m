function [matches] = f_track_neightopo_match2( part_A, part_B, f_o_s, n_neighbors)
%FUNCTION matches = f_track_neightopo_match2(part_A,part_B,f_o_s,n_neighbors)
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

%   INPUT:      part_A        - coordinates of particles in image A [n x 2]
%    
%               part_B        - coordinates of particles in image B [n x 2]
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
  

%%
% Calculate neighbourhood indices for each particle in image A and image B
for tempAorB = 1:2 
    
    if tempAorB==1 % Image A
        neighborInd = knnsearch(part_A(:,1:2),part_A(:,1:2),'K',n_neighbors+1);
        part_curr = part_A; 
    elseif tempAorB==2 % Image B
        neighborInd = knnsearch(part_B(:,1:2),part_B(:,1:2),'K',n_neighbors+1);
        part_curr = part_B;
    end
    
    featureVec_r = zeros(size(part_curr,1),n_neighbors);
    featureVec_phi = zeros(size(part_curr,1),n_neighbors);

    for parInd = 1:size(part_curr,1) % For each particle in part_curr

        % Grab neighboring particles of particle "parInd"
        parNeighborIndCoord = part_curr(neighborInd(parInd,:),:);

        % Compute distance btw neighboring particles to particle "parInd"
        dist_xy = parNeighborIndCoord(2:end,:) - repmat(part_curr(parInd,:), n_neighbors, 1);   
        r = sqrt( sum(dist_xy.^2,2) ) ; % abs distance 'r' is already in order 
        phi =  atan2( dist_xy(:,2) , dist_xy(:,1) ) ; % phi angle in xy-plane

        [ dist_phi_sort, dist_phi_sort_ind ] = sort(phi); % sort angle phi
        tempInd = find( dist_phi_sort_ind == 1 ); % find the index of r_min
        tempInd2 = [ tempInd : length(dist_phi_sort),  1:(tempInd-1) ]; % reorder index starting from the r_min

        r2 = r( dist_phi_sort_ind( tempInd2 ) ); % Reorder distance 'r' starting from r_min
        phi2 = phi( dist_phi_sort_ind( tempInd2 ) ); % Reorder angle 'phi' starting from r_min
        phiDiff2 = phi2([2:end,1]) - phi2; % Compute angle between two r's
        phiDiff2(phiDiff2 < 0) = phiDiff2(phiDiff2 < 0) + 2*pi; % Compensate 2*pi to make sure it's always positive

        featureVec_r(parInd,:) = [r2(:)']; % /min(r2); % Normalize r_min if needed
        featureVec_phi(parInd,:) = [phiDiff2(:)']; % Make sure the angle difference is a row array

    end
    
    % Store feature vectors
    if tempAorB==1 % Image A
        featureVec_r_A = featureVec_r; % feature vector in "r" dimension
        featureVec_phi_A = featureVec_phi; % feature vector in "phi" dim.,  sum of each row of featureVec_phi_A is 2*pi
    elseif tempAorB==2  % Image B
        featureVec_r_B = featureVec_r; % feature vector in "r" dimension
        featureVec_phi_B = featureVec_phi; % feature vector in "phi" dim., 
    end
    
end 


%% Find the most similar particle in part_B for each particle in part_A
waitbarOrNot = 0; % Waitbar to be used to check the code progres
if waitbarOrNot==1
    hbar = parfor_progressbar(size(part_A,1),'Computing...');  % create the progress bar in parallel computing "parfor"
end
clear matches;
if f_o_s == Inf
    neighborInd_InBNearA = 1:size(part_B,1); % Search the whole field
else 
    neighborInd_InBNearA_All = knnsearch(part_B(:,1:2),part_A(:,1:2),'K', n_neighbors); % Find pts in B near pts in A
end   
    
for parInd = 1:size(part_A,1) % For each particle in part_A
     
    if waitbarOrNot==1
        hbar.iterate(1);
    end
    
    %%%%% Only use neighbors within f_o_s %%%%%
    if f_o_s < Inf
        neighborInd_InBNearA = neighborInd_InBNearA_All(parInd,:);
        dist_neigh_InBNearA = sqrt(sum( (repmat(part_A(parInd,:),length(neighborInd_InBNearA),1) - part_B(neighborInd_InBNearA,:)).^2, 2));
        neighborInd_InBNearA = neighborInd_InBNearA(dist_neigh_InBNearA < sqrt(2)*f_o_s);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(neighborInd_InBNearA) == 0

        feaVec_r_A1 = featureVec_r_A(parInd,:); % figure, plot(feaVec_r_A1,'x');
        feaVec_phi_A1 = featureVec_phi_A(parInd,:);
        
        feaVec_r_B = featureVec_r_B(neighborInd_InBNearA(1:length(neighborInd_InBNearA)),:);
        feaVec_phi_B = featureVec_phi_B(neighborInd_InBNearA(1:length(neighborInd_InBNearA)),:);
          
        %%%%% LSQ difference between each topology feature %%%%%
        corr_r  = sum( (repmat(feaVec_r_A1,length(neighborInd_InBNearA),1)-feaVec_r_B).^2 , 2);
        corr_phi = sum( (repmat(feaVec_phi_A1,length(neighborInd_InBNearA),1)-feaVec_phi_B).^2 , 2);
         
        [~,ind_min_corr_r] = sort(corr_r);
        [~,ind_min_corr_phi] = sort(corr_phi);
          
        %%%%% Both "r-DIM" and "phi-DIM" features reach minimum %%%%%
        if (ind_min_corr_r(1)==ind_min_corr_phi(1)) % ||  (ind_min_corr_r(1)==ind_min_corr_phi(2))  ...
                 
            % discard displacements should be smaller than the field of search
            if sqrt( sum( ( part_A(parInd,:) - part_B(neighborInd_InBNearA(ind_min_corr_r(1)),:) ).^2 ) ) < f_o_s   
                matches{parInd} = [parInd, neighborInd_InBNearA(ind_min_corr_r(1))]; % save matches
            end
            
            %%%%% Here I comment sometimes there can be a relaxed criterion %%%%%  
            % elseif (ind_min_corr_r(2)==ind_min_corr_phi(1)) %%%JY!!!  || (ind_min_corr_r(2) == ind_min_corr_phi(2))
            % 
            %     % discard displacements should be smaller than the field of search
            %     if sqrt( sum( ( part_A(parInd,:) - part_B(neighborInd_InBNearA(ind_min_corr_r(2)),:) ).^2 ) ) < 0.5*f_o_s
            %         matches{parInd} = [parInd, neighborInd_InBNearA(ind_min_corr_r(2))]; % save matches
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

   