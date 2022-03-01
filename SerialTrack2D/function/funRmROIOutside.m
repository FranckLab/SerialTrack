function [uGrid,vGrid] = funRmROIOutside(xGrid,yGrid,imgMask,uGrid,vGrid)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply ImgRefMask to make u,v nans if there is a hole
imgSize = size(imgMask);
xGrid_temp = xGrid(:);
yGrid_temp = yGrid(:);
[row1,~] = find(xGrid_temp>1); [row3,~] = find(xGrid_temp<imgSize(1)-1);
[row2,~] = find(yGrid_temp>1); [row4,~] = find(yGrid_temp<imgSize(2)-1);
row = intersect(intersect(intersect(row1,row2),row3), row4);
rowOutside = setdiff([1:1:numel(xGrid)], row);

x0y0Ind = sub2ind(imgSize, round(xGrid(row)), round(yGrid(row)));
temp1 = imgMask(x0y0Ind);
temp1(isnan(temp1)) = 1;

uGrid(row(~temp1)) = nan;
vGrid(row(~temp1)) = nan;

uGrid(rowOutside) = nan;
vGrid(rowOutside) = nan;


