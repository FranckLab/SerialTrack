function [ matches ] = f_track_hist_match( part_A, part_B, f_o_s, n_neighbours, gauss_interp )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histogram matching tracking method based on:
% T. Fuchs, C.J. Kähler: "Volumetrische Messung wandnaher Strömungen"
% Fachtagung "Experimentelle Strömungsmechanik"
% 5.-7. September 2017, Karlsruhe
%
% This file is a copy from:
% T. Janke, R. Schwarze, K. Bauer: "Part2Track: A MATLAB package for double 
% frame and time resolved Particle Tracking Velocimetry". vol 11, 100413, 
% SoftwareX (2020). Link: https://github.com/Part2Track/Part2Track
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input:      part_A        - coordinates of particles in image A [n x 4]
%   ------
%               part_B        - coordinates of particles in image B [n x 4]
% 
%               f_o_s         - field of search [px]
%
%               n_neighbours  - number of neighbouring particles [integer]
%
%               gauss_interp  - flag to control refinement of displacement
%                               guess via gauss fitting of histogram
%                                  |
%                                  |--> = 0: no gauss fitting
%                                  |--> = 1: perform gauss fitting
%
%   Output:     matches     - list of indices of matching particles [m x 2]
%   -------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% origin: Thomas Janke / 21.09.17
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculate neighbourhood indices for each particle in image A
indice_neigh = knnsearch(part_A(:,1:2),part_A(:,1:2),'K',n_neighbours+1);

parfor ii=1:size(part_A,1)
%     disp(['Tracking of particle ',num2str(ii),'.']);
    
    % grab neigbouring particles
     part_neigbours = part_A(indice_neigh(ii,:),:);
    
	% calculate all possible displacements within field of search 
	% initialise distance variables
    dist_x=NaN; 
    dist_y=NaN;

    for jj=1:size(part_neigbours,1) % loop over all neighbours
        part_B_temp = part_B((abs(part_neigbours(jj,1)-part_B(:,1))< f_o_s) &...
                             (abs(part_neigbours(jj,2)-part_B(:,2))< f_o_s),: );

        dist_x = [dist_x;part_neigbours(jj,1)-part_B_temp(:,1)]; % collect distances in x-direction
        dist_y = [dist_y;part_neigbours(jj,2)-part_B_temp(:,2)]; % collect distances in y-direction                
    end

    
    bin_count = ceil(size(dist_x,1)/10); % determine bin numbers for histogram   
    if bin_count < 5
        bin_count = 5;      
    end
    
    [counts,centers] = hist(dist_x,bin_count); % create histogram of displacements in x

    [~, index_max] = max(counts); % search peak
    
    if index_max > 1 && index_max < bin_count
        switch gauss_interp % gauss interpolation
            case 1
                counts_log=log(counts(index_max-1:index_max+1)); % ln of counts
                delta_x = (counts_log(1)-counts_log(3))/(2*(counts_log(3)-2*counts_log(2)+counts_log(1))); % Gaussion interpolation
                dist_x_max = centers(index_max)+delta_x*diff(centers(1:2)); % Add sub-sample shift
            otherwise
                dist_x_max = centers(index_max); % x displacement with max. probability    
        end

        bin_count = ceil(size(dist_y,1)/10); % determine bin numbers for histogram
        if bin_count < 5
            bin_count = 5;
        end
        
        [counts,centers] = hist(dist_y,bin_count); % create histogram of displacements in y
        [~, index_max] = max(counts); % search peak     
        
        if index_max > 1 && index_max < bin_count
            switch gauss_interp % gauss interpolation
                case 1                     
                    counts_log=log(counts(index_max-1:index_max+1)); % ln of counts
                    delta_y = (counts_log(1)-counts_log(3))/(2*(counts_log(3)-2*counts_log(2)+counts_log(1))); % Gaussion interpolation
                    dist_y_max = centers(index_max)+delta_y*diff(centers(1:2)); % Add sub-sample shift
                otherwise
                    dist_y_max = centers(index_max); % x displacement with max. probability 
            end
                        
            temp_dist_x = part_A(ii,1)-part_B(:,1); % calculate x displacement of current particle to all particles in second frame
            temp_dist_y = part_A(ii,2)-part_B(:,2); % calculate y displacement of current particle to all particles in second frame 
            
            temp_dist_abs = sqrt(temp_dist_x.^2+temp_dist_y.^2); % calculate absolute displacement of current particle to all particles in second frame

            tempi = sqrt((temp_dist_x-dist_x_max).^2+(temp_dist_y-dist_y_max).^2); % calculate difference between most probable displacement and possible displacements
            [~, index_min] = min(tempi); % find displacement with minimum deviation 
            
            if min(tempi) < 5 % valid matches if difference between predicted and found displacements are smaller
                if temp_dist_abs(index_min) < f_o_s % discard displacements larger field of search
                    matches{ii} = [ii, index_min]; % save matches
                end
            end
        end
    end
end

% Check if there are found matches
if exist('matches','var')
    matches = cell2mat(matches');
else
    matches = [];
end

