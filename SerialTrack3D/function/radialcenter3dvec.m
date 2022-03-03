function [pts] = radialcenter3dvec(img, x0, beadParameter)
% [pts] = radialcenter3dvec(img, x0, beadParameter) localizes particles to
% sub-voxel accuracy using radial symmetry method based on
%
% Liu, Shu-Lin, et al. "Fast and high-accuracy localization for 
% three-dimensional single-particle tracking." Scientific reports 3 
% (2013): 2462.
%
% Note: Code modified from the radial symmetry code for 3D images provided 
% by Shulin Liu
%
% INPUTS
% -------------------------------------------------------------------------
%   u: cell containing the input displacement field. (u{1:3} = {u_x, u_y,
%   	u_z})
%   thr: theshold for passing residiual (default = 2)
%   epsilon: fluctuation level due to cross-correlation (default = 0.1)
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u: cell containing the displacement field with outliers removed
%   normFluctValues: normalized fluctuation values based on the universal
%   outier test.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Parse Inputs
winSize = beadParameter.winSize;
dccd = beadParameter.dccd;
dxccd = dccd(1); dyccd = dccd(2); dzccd = dccd(3);
abc = beadParameter.abc;
a = abc(1); b = abc(2); c = abc(3);
forloop = beadParameter.forloop;
randNoise = beadParameter.randNoise;

img = img + randNoise*rand(size(img));

%% Collect Neigbhorhood for each pt

% Pad Image with random numbers
nd = floor(winSize/2);
% img = padarray(img,nd);
img = padNoise(img,nd);

% Collecting neighbors for each pt
xtemp = bsxfun(@plus,x0,nd);    %Position change due to padding
I = zeros(winSize(1),winSize(2),winSize(3),size(x0,1));
for i = 1:size(x0,1)
    I(:,:,:,i) = img(xtemp(i,1)-nd(1):xtemp(i,1)+nd(1),xtemp(i,2)-nd(2):xtemp(i,2)+nd(2),...
        xtemp(i,3)-nd(3):xtemp(i,3)+nd(3));
end

%% Begin Localisation
[ny,nx,nz,nI] = size(I);
[px,py,pz] = meshgrid(-(nx-1)/2*dxccd:dxccd:(nx-1)/2*dxccd,-(ny-1)/2*dyccd:dyccd:(ny-1)/2*dyccd,-(nz-1)/2*dzccd:dzccd:(nz-1)/2*dzccd);

x=px(2:ny-1,2:nx-1,2:nz-1)/a;
y=py(2:ny-1,2:nx-1,2:nz-1)/b;
z=pz(2:ny-1,2:nx-1,2:nz-1)/c;

dx=dxccd/a;
dy=dyccd/b;
dz=dzccd/c;

xx = bsxfun(@times,I,px);
yy = bsxfun(@times,I,py);
zz = bsxfun(@times,I,pz);

bb = sum(sum(sum(I,1),2),3);
xm = sum(sum(sum(xx,1),2),3)./bb;
ym = sum(sum(sum(yy,1),2),3)./bb;
zm = sum(sum(sum(zz,1),2),3)./bb;

%% Calculate the gradients in x,y,z directions
d = sqrt((bsxfun(@minus,x,xm)).^2 + (bsxfun(@minus,y,ym)).^2 + (bsxfun(@minus,z,zm)).^2);

u = (I(2:ny-1,3:nx,2:nz-1,:) - I(2:ny-1,1:nx-2,2:nz-1,:))/dx;
v = (I(3:ny,2:nx-1,2:nz-1,:) - I(1:ny-2,2:nx-1,2:nz-1,:))/dy;
w = (I(2:ny-1,2:nx-1,3:nz,:) - I(2:ny-1,2:nx-1,1:nz-2,:))/dz;

dimg = sqrt((u.^2+v.^2+w.^2));

% weighting
q=dimg.^2./d;
u=u./dimg;
v=v./dimg;
w=w./dimg;

%% Calculate the center of the particle
A = zeros(3,3,nI);
B = zeros(3,nI);

%Get A matrix
A(1,1,:) = squeeze(sum(sum(sum(q.*(1-u.*u),1),2),3));
A(2,2,:) = squeeze(sum(sum(sum(q.*(1-v.*v),1),2),3));
A(3,3,:) = squeeze(sum(sum(sum(q.*(1-w.*w),1),2),3));
A(1,2,:) = -squeeze(sum(sum(sum(q.*u.*v,1),2),3));
A(2,1,:) = A(1,2,:);
A(1,3,:) = -squeeze(sum(sum(sum(q.*u.*w,1),2),3));
A(3,1,:) = A(1,3,:);
A(2,3,:) = -squeeze(sum(sum(sum(q.*v.*w,1),2),3));
A(3,2,:) = A(2,3,:);

%%get the B matrix
B(1,:) = -squeeze(sum(sum(sum(bsxfun(@times,q.*(u.*u-1),x) + ...
    bsxfun(@times,q.*u.*v,y) + bsxfun(@times,q.*u.*w,z),1),2),3));
B(2,:) = -squeeze(sum(sum(sum(bsxfun(@times,q.*u.*v,x) + ...
    bsxfun(@times,q.*(v.*v-1),y) + bsxfun(@times,q.*v.*w,z),1),2),3));
B(3,:) = -squeeze(sum(sum(sum(bsxfun(@times,q.*u.*w,x) + ...
    bsxfun(@times,q.*v.*w,y) + bsxfun(@times,q.*(w.*w-1),z),1),2),3));


if forloop == 0;
    
    % Create Sparse A matrix as Asparse
    t(1,1,:) = 3*(0:nI-1);
    AidxM = [1,2,3;1,2,3;1,2,3];
    AidxM = bsxfun(@plus,repmat(AidxM,[1,1,nI]),t);
    AidxN = [1,2,3;1,2,3;1,2,3]';
    AidxN = bsxfun(@plus,repmat(AidxN,[1,1,nI]),t);
    Asparse = sparse(AidxM(:),AidxN(:),A(:));
    
    B = B(:);
    
    %Solve
    xyz = Asparse\B;
    
elseif forloop == 1;
    
    %Solve
    xyz = zeros(nI*3,1);
    for i =1:nI
        xyz(3*i-2:3*i) = A(:,:,i)\B(:,i);
    end
    t(1,1,:) = 3*(0:nI-1);
end


% Switch to m,n,o format
pts(:,2) = xyz(t+1)*a;
pts(:,1) = xyz(t+2)*b;
pts(:,3) = xyz(t+3)*c;

nonNaNPts = sum(isnan(pts),2)<1;
pts = pts(nonNaNPts,:)+x0(nonNaNPts,:);

end


