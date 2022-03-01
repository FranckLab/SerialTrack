function [xGrid,yGrid,fGrid] = funScatter2Grid2D(varargin)
%FUNSCATTER2GRID2D: Interpolate 2D scattered data to gridded data.
%   [xGrid,yGrid,fGrid] = funScatter2Grid2D(x,y,f,sxy,smoothness)
%   
%   INPUT: 2D scatterred data with coordinates x,y and value f
%          sxy = [sx,sy] is the step for griddata
%          smoothness is the coefficient for regularization (function regularizeNd[1])
%          (Optional) gridded coordinates (xGrid,yGrid)
%
%   OUTPUT: 2D gridded data (fGrid) 
%           Gridded coordinates (xGrid,yGrid)
%
% -----------------------------------------------
% Author: Jin Yang (jyang526@wisc.edu)
% Date: 06-24-2020
%
% Reference
% [1] https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x,y,f,sxyz,smoothness,xGrid,yGrid] = parseargs(varargin);
x=x(:); y=y(:); smoothness=abs(smoothness);

if isempty(smoothness)
    alpha = 0; % No regularization
end
if isempty(sxyz)
    sxyz = floor((max([x,y])-min([x,y]))/20);
end
    
if isempty(xGrid) || isempty(yGrid)
    xList = min(x(:)):sxyz(1):max(x(:)+sxyz(1));
    yList = min(y(:)):sxyz(2):max(y(:)+sxyz(2));
    [yGrid,xGrid] = meshgrid(yList,xList);
else
    xList = unique(xGrid);  xList=xList(:)';
    yList = unique(yGrid);  yList=yList(:)';
end
 

% ------ Regularization ------
if smoothness==0
    f_interp = scatteredInterpolant(x,y,f,'linear');
    fGrid = reshape(f_interp(xGrid(:),yGrid(:)), size(xGrid));
else
    fGrid = regularizeNd([x,y],f,{xList,yList},smoothness);
end

% ------ Plot and check ------
% figure, scatter3(xGrid(:),yGrid(:),zGrid(:),ones(length(xGrid(:)),1),fGrid(:)); cb=colorbar;
% figure, scatter3(xGrid(:),yGrid(:),zGrid(:),ones(length(xGrid(:)),1),fGrid2(:)); cb=colorbar;
% figure, scatter3(x,y,z,ones(length(x(:)),1),f(:)); cb=colorbar;


end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [x,y,f,sxy,smoothness,xGrid,yGrid] = parseargs(vargin)
x=[]; y=[]; f=[]; sxy=[]; smoothness=[]; xGrid=[]; yGrid=[];
x=vargin{1}; y=vargin{2}; f=vargin{3};
try 
    sxy = vargin{4};
    try
        smoothness = vargin{5};
    catch
    end
catch
end
 
    
try
    xGrid=vargin{6};
    yGrid = vargin{7};
catch
end

end