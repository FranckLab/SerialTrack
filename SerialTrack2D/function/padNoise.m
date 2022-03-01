function [I] = padNoise(I,padSize)
I = padarray(I,padSize);

pds1 = 1:padSize(1); pds1 = [pds1,size(I,1)-padSize(1)+1:size(I,1)];
pds2 = 1:padSize(2); pds2 = [pds2,size(I,2)-padSize(2)+1:size(I,2)];
pds3 = 1:padSize(3); pds3 = [pds3,size(I,3)-padSize(3)+1:size(I,3)];

I0 = zeros(size(I));

I0(pds1,:,:) = 1;
I0(:,pds2,:) = 1;
I0(:,:,pds3) = 1;
idx = find(I0==1);

I(idx) = 0.00005*rand(size(idx)) + 0.00005;

end

