function [pts] = radial2center(I,x0,beadParameter)

winSize = beadParameter.winSize;  %Window size

%% Pad image with noise to take of partial beads at the boundary of image


nd = floor(winSize/2);

I = padarray(I,nd);
padSize = nd;
pds1 = 1:padSize(1); pds1 = [pds1,size(I,1)-padSize(1)+1:size(I,1)];
pds2 = 1:padSize(2); pds2 = [pds2,size(I,2)-padSize(2)+1:size(I,2)];

I0 = zeros(size(I));

I0(pds1,:) = 1;
I0(:,pds2) = 1;
idx = find(I0==1);
I(idx) = 0.0000005*rand(size(idx)) + 0.0000005;

%% Select regions around each particle

xtemp = bsxfun(@plus,x0,nd);    %Position change due to padding
I0 = zeros(winSize(1),winSize(2),size(x0,1));

for i = 1:size(x0,1)
    I0(:,:,i) = I(xtemp(i,1)-nd(1):xtemp(i,1)+nd(1),xtemp(i,2)-nd(2):xtemp(i,2)+nd(2));
end

pts = zeros(size(x0));
for i = 1:size(x0,1)
    [pts(i,2),pts(i,1),~] = radialcenter(I0(:,:,i));
end

nonNaNPts = sum(isnan(pts),2)<1;
pts = pts(nonNaNPts,:)+x0(nonNaNPts,:);

% pts = bsxfun(@plus,pts,-nd+1); % Arrange for center position shift from padding

end

