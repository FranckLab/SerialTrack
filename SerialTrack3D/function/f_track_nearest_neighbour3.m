function [ matches ] = f_track_nearest_neighbour3( part_A, part_B, f_o_s )
%FUNCTION matches = f_track_nearest_neighbour3(part_A, part_B, f_o_s)
% Objective: Nearest neigbhour tracking routine.
% ---------------------------------------------------
%
%   INPUT:      part_A      - coordinates of particles in image A [n x 3]
%    
%               part_B      - coordinates of particles in image B [n x 3]
% 
%               f_o_s       - field of search [px]
%
%   OUTPUT:     matches     - list of indices of matching particles [m x 2]
%    
% ---------------------------------------------------
% References
% [1] T Janke, R Schwarze, K Bauer. Part2Track: A MATLAB package for double
%     frame and time resolved Particle Tracking Velocimetry. 11, 100413, SoftwareX (2020).
% ---------------------------------------------------
% origin: Thomas Janke / 2017.09
% modified: Jin Yang / 2020.12
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
for parInd = 1:size(part_A,1) % Loop all particles in Image A
      
    temp_dist = sqrt((part_B(:,1)-part_A(parInd,1)).^2 + (part_B(:,2)-part_A(parInd,2)).^2 + ...
        (part_B(:,3)-part_A(parInd,3)).^2 ); % Calculate all possible distances between particle parInd and particles in image B
    
    [~,index_min] = min(temp_dist); % Find candidate with minimum displacement
    matches{parInd} = [parInd, index_min];
         
end

% Check if there are found matches
if exist('matches','var')
    matches = cell2mat(matches');
else
    matches = [];
end


% parMatchInd = 1; % Index of matched particle-pair  
% hbar = waitbar(0,'Please wait.');
% for parInd = 1:size(part_A,1) % Loop all particles in Image A
%      
%     waitbar(parInd/size(part_A,1));
%     
%     temp_dist = sqrt((part_B(:,1)-part_A(parInd,1)).^2 + (part_B(:,2)-part_A(parInd,2)).^2 + ...
%         (part_B(:,3)-part_A(parInd,3)).^2 ); % Calculate all possible distances between particle parInd and particles in image B
%     
%     if min(temp_dist) < f_o_s % Just consider candidates within field of search
%         
%         [~,index_min] = min(temp_dist); % Find candidate with minimum displacement
%         matches(parMatchInd,1:2) = [parInd, index_min];
%         
%         parMatchInd = parMatchInd+1;
%         
%     end
% end
% 
% close(hbar);

end

