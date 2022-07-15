function [part_list] = f_detect_particles3(im,beadPara)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Particle detection using Laplacian of Gaussian image filtering based on:
% Joris Heyman TracTrac: "A fast multi-object tracking algorithm for motion
% estimation", Computer & Geosciences, 2019, 128, 11-18
% Github: https://github.com/jorishey1234/tractrac
% under MIT License (see below)
% Copyright (c) 2018 jorishey1234
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input:      im            - image 
%   ------
%               p_size        - particle size [px]
% 
%               p_int         - particle intensity [counts]
%
%   Output:     part_list     - list of detected particles [m x 2]
%   -------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p_size = beadPara.beadRad;
p_int = beadPara.thres;
  
%% Find beads between minPixels and maxPixels

%Parameters
thres = beadPara.thres;    %Threshold value
minPixels = beadPara.minSize;  %Minimum pixel count in blob for bead
maxPixels = beadPara.maxSize;  %Maximum pixel count in blob for bead

%Binary Image % figure, imagesc3(im);
BW = im>thres; % figure, imshow(BW);

% Find bead blobs
CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
beadBlob = numPixels>minPixels & numPixels<maxPixels;

% Find beads between minPixels and maxPixels
BW1 = bwareaopen(BW, minPixels );  % figure, imshow(BW1);
BW2 = bwareaopen(BW, maxPixels );  % figure, imshow(BW2);
BW3 = BW1-BW2; % figure, imshow(BW3);
if sum(mean(BW3(:))) < 1, BW3 = BW1; end 
    
im = im.*BW3; % figure, imagesc3(im);
% figure, imagesc3(im); im = im.*BW3; 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if p_size > 0
 
    % Create image filter (Laplacian of Gaussian 'log')
    flog = fspecial3('log',[1 1 1]*ceil(p_size)*2+1, p_size) ;

    % Apply LoG filter
    im_filt=imfilter(im,-flog,'replicate','same'); % figure, imagesc3(im_filt);

    %% Find particle centers within pixel accuracy

    r=rand(size(im_filt))*1e-5; % add small random noise
    SE = strel('sphere',(p_size*2)+1); % create square element larger than p_size

    % local maximum filter
    b = imdilate(im_filt+r,SE); % apply SE on image

    % find maxima above threshold
    k = find((im_filt+r==b)&(im>p_int));
    [y,x,z] = ind2sub(size(im_filt),k);
    
    % Remove points on border
    w=size(im_filt,2); h=size(im_filt,1); t=size(im_filt,3);
    nb=floor((p_size+2)/2);
    id=find(y<h-nb & y>=nb+1 & x<w-nb & x>=nb+1 & z<t-nb & z>=nb+1);
    x=x(id); y=y(id); z=z(id);

    %% Subpixel Interpolation
    % Gaussian interpolation
    % dx,dy - subpixel shifts in x & y direction
    [dx,dy,dz,~]=subpix2nd3(real(log(im_filt-min(im_filt(:))+1e-8)),x,y,z,p_size);

    % Take only peaks that moved into the pixel
    indice=find( abs(dx)<0.5 & abs(dy)<0.5 & abs(dz)<0.5 );
    xb=x+dx; % add sub-pixel position in x
    yb=y+dy; % add sub-pixel position in y
    zb=z+dz; % add sub-pixel position in z

    xm=xb(indice);
    ym=yb(indice);
    zm=zb(indice);

    % create final output matrix
    part_list = [ym xm zm];
    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif p_size == 0
    
    
    %%%%% regionprops %%%%%
    CC = bwconncomp(BW3);
    props = regionprops3( CC ,'Volume','Centroid','VoxelIdxList','PrincipalAxisLength');

    propsVolList = [props.Volume]';
    propsRadList = sort( (propsVolList/4*3/pi).^(1/3), 'descend');
    propMajorAxisLength = [props.PrincipalAxisLength(:,1)]';
    propMinorAxisLength = [props.PrincipalAxisLength(:,3)]';
    propCentroidList = [props.Centroid]'; 
     
    propsEccentricity = propMajorAxisLength ./ propMinorAxisLength;

    % figure, hist(propsRadList,20);
    % figure, hist(propsEccentricity,20);
    
    [row1,col1] = find( propsRadList<maxPixels+2 & propsRadList>minPixels-1 & ...
                        propsEccentricity > 0.25 & propsEccentricity < 4);
    
    part_list = propCentroidList([2,1,3],col1)';
    % figure, imshow(BW3); hold on; plot(propCentroidList(:,1),propCentroidList(:,2),'bo');
    

end
    
end




%% 
function [Dx,Dy,Dz,Z]=subpix2nd3(a,x,y,z,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Particle detection within sub-pixel accuracy based on:
% Joris Heyman TracTrac: "A fast multi-object tracking algorithm for motion
% estimation", Computer & Geosciences, 2019, 128, 11-18
% from: https://github.com/jorishey1234/tractrac
% under MIT License (see below)
% Copyright (c) 2018 jorishey1234
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

np=floor(n/2); % number of points for pencil
pencil=-np:np; % create pencil
X=repmat(pencil,length(x),1); % create pencil for every particle

n=length(pencil); % pencil size
YV=[]; % Vertical
YH=[]; % Horizontal
YZ=[]; % Thickness

% Convert indices
for i=1:length(pencil) 
idV = sub2ind(size(a),max(1,min(size(a,1),y+pencil(i))),x,max(1,min(size(a,3),z+pencil(i))));
idH = sub2ind(size(a),y,max(1,min(size(a,2),x+pencil(i))),max(1,min(size(a,3),z+pencil(i))));
idZ = sub2ind(size(a),max(1,min(size(a,1),y+pencil(i))),max(1,min(size(a,2),x+pencil(i))),z);
YV(:,i)=a(idV);
YH(:,i)=a(idH);
YZ(:,i)=a(idZ);
end

% 2nd order poly a+bx+cx^2=0
s2=sum(pencil.^2);
s4=sum(pencil.^4);

bH=sum((YH.*X),2)'./s2;
cH=-(-s2*sum(YH,2)'+n*sum(X.^2.*YH,2)')/(s2^2-s4*n);

bV=sum(YV.*X,2)'./s2;
cV=-(-s2*sum(YV,2)'+n*sum(X.^2.*YV,2)')/(s2^2-s4*n);

bZ=sum(YZ.*X,2)'./s2;
cZ=-(-s2*sum(YZ,2)'+n*sum(X.^2.*YZ,2)')/(s2^2-s4*n);

% Peaks on hor, vert, and z-axis
dH=-bH./2./cH;
dV=-bV./2./cV;
dZ=-bZ./2./cZ;
    
Dx=(dH');
Dy=(dV');
Dz=(dZ');
Z=YH; %Not sure what it is
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIT License
% Copyright (c) 2018 jorishey1234
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

