function [XYZ,U,F] = funCompDefGrad3(uvw, part_A, f_o_s, n_neighbors)
%FUNCTION [XYZ,U,F] = funCompDefGrad(uvw, part_A, f_o_s, n_neighbors)
% Objective: To compute tracked displacements and def grad tensor based on
%            moving least square fitting
% ---------------------------------------------------
% 
%   INPUT:  uvw                     Tracked displacement components [n x 3]
%           parCoordA               Coordinates of particles in image A
%           f_o_s                   field of search [px]
%           n_neighbours            number of neighbouring particles [integer]
 
%   OUTPUT: XYZ                     Coordinates of particles in image A
%           U                       Disp vectoc U = [U_node1, V_node1, W_node1, ..., U_nodeN, V_nodeN, W_nodeN]';
%           F                       Deformation gradient tensor is assembled into a long vector:
%                                   F = [F11_node1, F21_node1, F31_node1, ...
%                                        F12_node1, F22_node1, F32_node1, ... 
%                                        F13_node1, F23_node1, F33_node1, ... 
%                                        ..., 
%                                        F11_nodeN, F21_nodeN, F31_nodeN, ... 
%                                        F12_nodeN, F22_nodeN, F32_nodeN, ... 
%                                        F13_nodeN, F23_nodeN, F33_nodeN]';
% 
% ---------------------------------------------------
% Author: Jin Yang
% Contact and support: jyang526@wisc.edu
% Date: 2020.12.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Example:
% uvw = uAB_; % displacement: [u,v,w]
% part_A = partA(matchesAB(:,1),:);
% part_B = partB(matchesAB(:,2),:);
% n_neighbors = 20;
% f_o_s = 50;
 

%% Initialization
U = nan(size(uvw,1)*3,1); F = repmat(U,3,1); XYZ = nan(size(uvw,1),3);

%%
% Calculate neighbourhood indices for each particle in image A
neighborInd = knnsearch(part_A(:,1:3),part_A(:,1:3),'K',n_neighbors+1);

for parInd = 1:size(part_A,1)
    
    parNeighborInd = neighborInd(parInd,:);
    parNeighborIndCoord = part_A(parNeighborInd,:);
    
    try 
    
        % %%%%% Only use neighbors within f_o_s %%%%%
        % Compute distance
        dist_xyz = parNeighborIndCoord - repmat(part_A(parInd,:), n_neighbors+1, 1);
        r =  sqrt( sum(dist_xyz.^2,2) ) ; % Abs distance r: already in order
        
        parNeighborInd = parNeighborInd(r < sqrt(3)*f_o_s);
   
        % Coordinates: part_A(parNeighborInd, 1:3)
        % u: uvw(parNeighborInd, 1)
        % v: uvw(parNeighborInd, 2)
        % w: uvw(parNeighborInd, 3)
        % 
        % % We need at least four data points to compute F={u,x v,x w,x u,y,v,y w,y u,z v,z w,z}
        %   x = x0 + u + u,x*(x-x0) + u,y*(y-y0) + u,z*(z-z0) 
        %   y = y0 + v + v,x*(x-x0) + v,y*(y-y0) + v,z*(z-z0) 
        %   z = z0 + w + w,x*(x-x0) + w,y*(y-y0) + w,z*(z-z0) 
        
        AMatrix = [ones(length(parNeighborInd),1), part_A(parNeighborInd,1)-part_A(parInd,1), ...
             part_A(parNeighborInd,2)-part_A(parInd,2), part_A(parNeighborInd,3)-part_A(parInd,3)];
         
        UPara = AMatrix \ uvw(parNeighborInd,1);
        VPara = AMatrix \ uvw(parNeighborInd,2);
        WPara = AMatrix \ uvw(parNeighborInd,3);

        U(3*parInd-2:3*parInd) = [UPara(1);VPara(1);WPara(1)];
        F(9*parInd-8:9*parInd) = reshape([UPara(2:4)'; VPara(2:4)'; WPara(2:4)'],9,1); 
        XYZ(parInd,1:3) = part_A(parInd,1:3);
        
    catch
 
        
        U(3*parInd-2:3*parInd) = [uvw(parInd,1);uvw(parInd,2);uvw(parInd,3)];
        F(9*parInd-8:9*parInd) = nan(9,1);
        XYZ(parInd,1:3) = part_A(parInd,1:3);
        
    end
    
           
       
       

end