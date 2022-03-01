function [track] = removeOutlierTPT2(varargin)
%REMOVEOUTLIERTPT2: to find outliers in scattered data based on a residual method [1].
%  [track] = removeOutlierPSPT(varargin) 
%
%
% INPUTS
% -------------------------------------------------------------------------
%   x0:     Scattered point in nD scape at t=0;
%   x1:     Scattered point in nD scape at t=1;
%   track:  Track results
%   thres:  Threshold value for passing residual
%
% OUTPUTS
% -------------------------------------------------------------------------
%   track:  Only valid track results saved. Outliers are removed
%
%
% -------------------------------------------------------------------------
% Reference
% [1] Westerweel, Jerry, and Fulvio Scarano. "Universal outlier detection for 
% PIV data." Experiments in fluids 39.6 (2005): 1096-1100.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%% Parse Inputs %%%%%%%%%%%
x0 = varargin{1};
x1 = varargin{2};
track = varargin{3};
thres = varargin{4};

%Number of neighbors to compute median and residual displacement
nNeigh = 40; 


%%%%%%%%%%% Find Outliers %%%%%%%%%%%
u = x1(track(track>0),:) - x0(track>0,:);    %Displacement
u1 = u(:,1); u2 = u(:,2);  % Disp components

x = x0(track>0,:);      % Data points for which displacement tracked. 
idx = knnsearch(x,x,'k',nNeigh);
% idx = idx(:,2:end); % Remove first point. It's the self point
un1 = zeros(size(idx)); un2 = un1;  

%Find median u (um)
un1(:) = u1(idx(:)); um1 = median(un1,2);
un2(:) = u2(idx(:)); um2 = median(un2,2);
 

%Find residual
e = 0.1;
r1 = abs(bsxfun(@plus,un1,-um1)); r1 = median(r1,2)+e;
r2 = abs(bsxfun(@plus,un2,-um2)); r2 = median(r2,2)+e;
 

%Normalize residual
rn1 = abs(bsxfun(@plus,u1,-um1)); rn1 = rn1./r1;
rn2 = abs(bsxfun(@plus,u2,-um2)); rn2 = rn2./r2;
 

% Threshold out values
ridx1 = rn1>thres; 
ridx2 = rn2>thres;
 
IDX = ridx1 | ridx2;

% Remove wrong IDX from track
idx = find(track>0);
track(idx(IDX)) = 0;

end

