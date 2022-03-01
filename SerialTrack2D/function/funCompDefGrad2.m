function [XY,U,F] = funCompDefGrad2(uv, parCoord, f_o_s, n_neighbors)
%FUNCOMPDEFGRAD2: to compute tracked displacements and deformation gradient
%                 tensor based on the moving least square fitting method
%   [XY,U,F] = funCompDefGrad2(uvw, part_A, f_o_s, n_neighbors) 
% ---------------------------------------------------
% 
%   INPUT:  uv                      Tracked displacement components [n x 2]
%           parCoord                Coordinates of particles
%           f_o_s                   field of search [px]
%           n_neighbours            number of neighbouring particles [integer]
 
%   OUTPUT: XY                      Coordinates of particles in image A
%           U                       Disp vectoc U = [U_node1, V_node1, ..., U_nodeN, V_nodeN]';
%           F                       Deformation gradient tensor is assembled into a long vector:
%                                   F = [F11_node1, F21_node1, F12_node1, F22_node1,  
%                                        ..., 
%                                        F11_nodeN, F21_nodeN, F12_nodeN, F22_nodeN]';
% 
% ---------------------------------------------------
% Author: Jin Yang
% Contact and support: jyang526@wisc.edu
% Date: 2020.12.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Example:
% uv = uAB_; % displacement: [u,v]
% part_A = partA(matchesAB(:,1),:);
% part_B = partB(matchesAB(:,2),:);
% n_neighbors = 20;
% f_o_s = 50;
 

%% Initialization
U = nan(size(uv,1)*2,1); F = repmat(U,2,1); XY = nan(size(uv,1),2);

%%
% Calculate neighbourhood indices for each particle in image A
neighborInd = knnsearch(parCoord(:,1:2),parCoord(:,1:2),'K',n_neighbors+1);

for parInd = 1:size(parCoord,1)
    
    parNeighborInd = neighborInd(parInd,:);
    parNeighborIndCoord = parCoord(parNeighborInd,:);
    
     try 
    
        % %%%%% Only use neighbors within f_o_s %%%%%
        % Compute distance
        dist_xy = parNeighborIndCoord - repmat(parCoord(parInd,:), size(parNeighborIndCoord,1), 1);
        r =  sqrt( sum(dist_xy.^2,2) ) ; % Abs distance r: already in order
        
        parNeighborInd = parNeighborInd(r < 2*f_o_s);
   
        % Coordinates: part_A(parNeighborInd, 1:2)
        % u: uvw(parNeighborInd, 1)
        % v: uvw(parNeighborInd, 2)
        %  
        % % We need at least four data points to compute F = {u,x v,x u,y v,y}
        %   x = x0 + u + u,x*(x-x0) + u,y*(y-y0)   
        %   y = y0 + v + v,x*(x-x0) + v,y*(y-y0)  
         
        AMatrix = [ones(length(parNeighborInd),1), parCoord(parNeighborInd,1)-parCoord(parInd,1), ...
             parCoord(parNeighborInd,2)-parCoord(parInd,2) ];
         
        UPara = AMatrix \ uv(parNeighborInd,1);
        VPara = AMatrix \ uv(parNeighborInd,2);
        
        U(2*parInd-1:2*parInd) = [UPara(1);VPara(1)];
        F(4*parInd-3:4*parInd) = reshape([UPara(2:3)'; VPara(2:3)'],4,1); 
        XY(parInd,1:2) = parCoord(parInd,1:2);
        
    catch
  
        U(2*parInd-1:2*parInd) = [uv(parInd,1);uv(parInd,2)];
        F(4*parInd-3:4*parInd) = nan(4,1);
        XY(parInd,1:2) = part_A(parInd,1:2);
        
    end
    
           
       
       

end



