function [coordinatesFEM,elementsFEM] = funMeshSetUp(x,y) 

M = size(x,1);  N = size(x,2);   % N is vertically in image; M is horizontally in image;
coordinatesFEM = zeros(M*N ,2);

 
% I have transpose x and y because Matlab is read matrix in column direction
for i = 1:size(coordinatesFEM,1)
    coordinatesFEM(i,:)  = [x(i),y(i)];
    % x is horizontal position in the image
    % y is vertical position in the image
end

elementsFEM = zeros((M-1)*(N-1),4);
for j = 1:N-1
    for i = 1:M-1
        elementsFEM((j-1)*(M-1)+i ,:) = [(j-1)*(M)+i (j-1)*(M)+i+1 j*(M)+i+1 j*(M)+i];
    end
end