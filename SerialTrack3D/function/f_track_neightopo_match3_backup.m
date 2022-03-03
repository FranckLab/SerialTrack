function [matches] = f_track_neightopo_match3_backup( part_A, part_B, f_o_s, n_neighbors)
%FUNCTION matches = f_track_neightopo_match3(part_A,part_B,f_o_s,n_neighbors)
% Objective: Tracking based on topology similarity
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
% Calculate neighbourhood indices for each particle in image A
neighborInd = knnsearch(part_A(:,1:3),part_A(:,1:3),'K',n_neighbors+1);
part_curr = part_A;
featureVec_r = zeros(size(part_curr,1),n_neighbors);
featureVec_phi = zeros(size(part_curr,1),n_neighbors);
featureVec_theta = zeros(size(part_curr,1),n_neighbors);
 
for parInd = 1:size(part_curr,1) % For each particle in part_A
   
    % Grab neighboring particles  
    parNeighborIndCoord = part_curr(neighborInd(parInd,:),:);
    
    % Compute distance 
    dist_xyz = parNeighborIndCoord(2:end,:) - repmat(part_curr(parInd,:), n_neighbors, 1);   
    r = sqrt( sum(dist_xyz.^2,2) ) ; % Abs distance r: already in order 
    phi =  atan2( dist_xyz(:,2) , dist_xyz(:,1) ) ; % phi angle in xy-plane
    theta = atan2( dist_xyz(:,3), sqrt( sum(dist_xyz(:,1:2).^2,2) )  ); % polar angle in z-axis
    
    [ dist_phi_sort, dist_phi_sort_ind ] = sort(phi);
    tempInd = find( dist_phi_sort_ind == 1 ); % Find the index of min dist_r
    tempInd2 = [ tempInd : length(dist_phi_sort),  1:(tempInd-1) ]; % Reorder index starting from min dist_r
    
    r2 = r( dist_phi_sort_ind( tempInd2 ) ); % Reordered dist
    phi2 = phi( dist_phi_sort_ind( tempInd2 ) );
    phiDiff2 = phi2([2:end,1]) - phi2;
    phiDiff2(phiDiff2 < 0) = phiDiff2(phiDiff2 < 0) + 2*pi; % Reordered angle difference in xy-plane
    theta2 = theta( dist_phi_sort_ind( tempInd2 )  ); % Reordered angle in z-axis
     
    featureVec_r(parInd,:) = [r2];
    featureVec_phi(parInd,:) = [phiDiff2];
    featureVec_theta(parInd,:) = [theta2];
     
end

featureVec_r_A = featureVec_r;
featureVec_phi_A = featureVec_phi;
featureVec_theta_A = featureVec_theta;


%%
% Calculate neighbourhood indices for each particle in image B
neighborInd = knnsearch(part_B(:,1:3),part_B(:,1:3),'K',n_neighbors+1);
part_curr = part_B;
featureVec_r = zeros(size(part_curr,1),n_neighbors);
featureVec_phi = zeros(size(part_curr,1),n_neighbors);
featureVec_theta = zeros(size(part_curr,1),n_neighbors);

for parInd = 1:size(part_curr,1) % For each particle in part_A
   
    % Grab neighboring particles  
    parNeighborIndCoord = part_curr(neighborInd(parInd,:),:);
    
    % Compute distance 
    dist_xyz = parNeighborIndCoord - repmat(part_curr(parInd,:), n_neighbors+1, 1);  dist_xyz=dist_xyz(2:end,:);
    r =  sqrt( sum(dist_xyz.^2,2) ) ; % Abs distance r: already in order 
    phi =  atan2( dist_xyz(:,2) , dist_xyz(:,1) ) ; % phi angle in xy-plane
    theta = atan2( dist_xyz(:,3), sqrt( sum(dist_xyz(:,1:2).^2,2) )  ); % polar angle in z-axis
    
    [ dist_phi_sort, dist_phi_sort_ind ] = sort(phi);
    tempInd = find( dist_phi_sort_ind == 1 ); % Find the index of min dist_r
    tempInd2 = [ tempInd : length(dist_phi_sort),  1:(tempInd-1) ]; % Reorder index starting from min dist_r
    
    r2 = r( dist_phi_sort_ind( tempInd2 ) ); % Reordered dist
    phi2 = phi( dist_phi_sort_ind( tempInd2 ) );
    phiDiff2 = phi2([2:end,1]) - phi2;
    phiDiff2(phiDiff2 < 0) = phiDiff2(phiDiff2 < 0) + 2*pi; % Reordered angle difference in xy-plane
    theta2 = theta( dist_phi_sort_ind( tempInd2 )  ); % Reordered angle in z-axis
     
    featureVec_r(parInd,:) = [r2(:)'];
    featureVec_phi(parInd,:) = [phiDiff2(:)'];
    featureVec_theta(parInd,:) = [theta2(:)'];
     
end

featureVec_r_B = featureVec_r;
featureVec_phi_B = featureVec_phi;
featureVec_theta_B = featureVec_theta;


%% Find the max similar particle in part_B for each particle in part_A
% hbar = parfor_progressbar(size(part_A,1),'Computing...');  %create the progress bar
clear matches
if f_o_s == Inf
    neighborInd_InBNearA = 1:size(part_B,1); % Search the whole field
else 
    neighborInd_InBNearA_All = knnsearch(part_B(:,1:3),part_A(:,1:3),'K', n_neighbors); % Find pts in B near pts in A
end   
    
for parInd = 1:size(part_A,1) % For each particle in part_A
     
    % hbar.iterate(1);
    % disp(['Tracking of particle ',num2str(parInd),'.']);
     
    % %%%%% Only use neighbors within f_o_s %%%%%
    if f_o_s < Inf
        neighborInd_InBNearA = neighborInd_InBNearA_All(parInd,:);
     
        dist_neigh_InBNearA = sqrt(sum( (repmat(part_A(parInd,:),length(neighborInd_InBNearA),1) - part_B(neighborInd_InBNearA,:)).^2, 2));
        neighborInd_InBNearA = neighborInd_InBNearA(dist_neigh_InBNearA < sqrt(3)*f_o_s);
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

        % %%%%% Two of the three features reach minimum %%%%% 
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
        % %%%%% All three features reach minimum %%%%%
        if (ind_min_corr_r==ind_min_corr_phi) && (ind_min_corr_phi==ind_min_corr_theta)

            
            % discard displacements should be smaller than the field of search
            if sqrt( sum( ( part_A(parInd,:) - part_B(neighborInd_InBNearA(ind_min_corr_r),:) ).^2 ) ) < 0.5*f_o_s   
                matches{parInd} = [parInd, neighborInd_InBNearA(ind_min_corr_r)]; % save matches
            end
        end
    
     end % END of if isempty(neighborInd_BA) == 0
        
      
end
    

% Check if there are found matches
if exist('matches','var')
    matches = cell2mat(matches');
else
    matches = [];
end

% close(hbar)

  

%% %%%%% Old codes to plot each feature: corr_r, corr_phi, and corr_theta %%%%%
% close all
% parInd = 1800;
% feaVec_r_A1 = featureVec_r_A(parInd,:);
% 
% neighborInd_AB = knnsearch(part_B(:,1:3),part_A(parInd,1:3),'K',n_neighbors); % Find pts in B near pts in A
% 
% feaVec_r_B = [];
% for tempi = 1:n_neighbors
% feaVec_r_B(tempi,:) = featureVec_r_B(neighborInd_AB(tempi),:);
% end
%  
% figure, plot(feaVec_r_A1,'ko-','linewidth',2); 
% for tempi = 1:n_neighbors
% hold on; plot(feaVec_r_B(tempi,:),'+-');
% end
% 
% corr = [];
% for tempi = 1:n_neighbors
% corr(tempi) = sum((feaVec_r_A1-feaVec_r_B(tempi,:)).^2); %max(xcorr(feaVec_r_A1, feaVec_r_B(tempi,:)));
% end
% 
% figure, plot(corr,'o')
% 
% 
% %%
% feaVec_phi_A1 = featureVec_phi_A(parInd,:);
% 
% neighborInd_AB = knnsearch(part_B(:,1:3),part_A(parInd,1:3),'K',n_neighbors); % Find pts in B near pts in A
% 
% feaVec_phi_B = [];
% for tempi = 1:n_neighbors
% feaVec_phi_B(tempi,:) = featureVec_phi_B(neighborInd_AB(tempi),:);
% end
%  
% figure, plot(feaVec_phi_A1,'ko-','linewidth',2); 
% for tempi = 1:n_neighbors
% hold on; plot(feaVec_phi_B(tempi,:),'+-');
% end
% 
% corr = [];
% for tempi = 1:n_neighbors
% corr(tempi) = sum((feaVec_phi_A1-feaVec_phi_B(tempi,:)).^2); %max(xcorr(feaVec_phi_A1, feaVec_phi_B(tempi,:)));
% end
% 
% figure, plot(corr,'o')
% 
% 
% %%
% feaVec_theta_A1 = featureVec_theta_A(parInd,:);
% 
% neighborInd_AB = knnsearch(part_B(:,1:3),part_A(parInd,1:3),'K',n_neighbors); % Find pts in B near pts in A
% 
% feaVec_theta_B = [];
% for tempi = 1:n_neighbors
% feaVec_theta_B(tempi,:) = featureVec_theta_B(neighborInd_AB(tempi),:);
% end
%  
% figure, plot(feaVec_theta_A1,'ko-','linewidth',2); 
% for tempi = 1:n_neighbors
% hold on; plot(feaVec_theta_B(tempi,:),'+-');
% end
% 
% corr = [];
% for tempi = 1:n_neighbors
% corr(tempi) = sum((feaVec_phi_A1-feaVec_phi_B(tempi,:)).^2); %max(xcorr(feaVec_theta_A1, feaVec_theta_B(tempi,:)));
% end
% 
% figure, plot(corr,'o')






