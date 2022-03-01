function [track_A2B, u_A2B] = funCompDisp2(parCoordA,parCoordB,matches_A2B,outlrThres)
%FUNCOMPDISP2: to compute tracked displacements from matched particle pairs
%              and remove outliers
% FUNCTION [track_A2B, u_A2B] = funCompDisp2(parCoordA, parCoordB, matches_A2B, removeOutlierOrNot) 
% ---------------------------------------------------
% 
%   INPUT: parCoordA                Coordinates of particles in image A
%          parCoordB                Coordinates of particles in image A
%          matches_A2B              Matched particle pairs between image A and image B
%          outlrThres               To remove outliers if outlrThres>0
%
%   OUTPUT: track_A2B               Indices of tracked particle index in image B for each particle (row#) in image A        
%           u_A2B                   Tracked displacement components [n x 2]
% 
%
%
% ---------------------------------------------------
% Author: Jin Yang
% Contact and support: jyang526@wisc.edu
% Date: 2020.12.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% ===== Compute displacement from tracked particles =====
parCoordA_ = parCoordA(matches_A2B(:,1),:); parCoordB_ = parCoordB(matches_A2B(:,2),:);
u_A2B = parCoordB_ - parCoordA_;  % displacement

%%%%% Plot %%%%%
% figure, quiver3(parCoordA_(:,1),parCoordA_(:,2),parCoordA_(:,3),u_A2B(:,1),u_A2B(:,2),u_A2B(:,3));
% set(gca,'fontsize',20); view(2); axis equal; axis tight; view(3);
 
%%%%% Assemble variable "track_A2B" %%%%%
track_A2B = zeros(size(parCoordA,1),1);
for tempi = 1:size(matches_A2B)
    track_A2B(matches_A2B(tempi,1)) = matches_A2B(tempi,2);  
end

%%%%% Remove outliers %%%%%
try
    if outlrThres > 0 % Remove Outliers
       [track_A2B] = removeOutlierTPT2(parCoordA,parCoordB,track_A2B,outlrThres);
    end
catch
end

% ===== Compute displacement from good tracked particles =====
[parAtrack,~] = find(track_A2B>0);
parCoordA_ = parCoordA(parAtrack,:); parCoordB_ = parCoordB(track_A2B(parAtrack),:);
u_A2B = parCoordB_ - parCoordA_;  % displacement

% figure, quiver3(parCoordA_(:,1),parCoordA_(:,2),parCoordA_(:,3),u_A2B(:,1),u_A2B(:,2),u_A2B(:,3));
% set(gca,'fontsize',20); view(2); axis equal; axis tight; view(3);


end