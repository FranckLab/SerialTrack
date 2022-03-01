function [NeumannBCInd_F, NeumannBCInd_u]  = funFDNeumannBCInd(size1coordinatesFEM,M,N,Rad)
% Find indices of finite difference scheme Neumann boundary points 

%%
temp = 1:1:size1coordinatesFEM;  

temp = temp'; temp = reshape(temp,M,N); 
temp2 = temp(Rad+1:M-Rad, Rad+1:N-Rad); temp2 = reshape(temp2, (M-2*Rad)*(N-2*Rad),1);
temp3 = zeros(4*(M-2*Rad)*(N-2*Rad),1);
for i = 1:(M-2*Rad)*(N-2*Rad), temp3(4*i-3:4*i) = [4*temp2(i)-3; 4*temp2(i)-2; 4*temp2(i)-1; 4*temp2(i)]; end
%atemp = a(temp3);
temp4 = zeros(2*(M-2*Rad)*(N-2*Rad),1);
for i = 1:(M-2*Rad)*(N-2*Rad), temp4(2*i-1:2*i) = [2*temp2(i)-1; 2*temp2(i)]; end
%btemp = b(temp4);

NeumannBCInd_F = temp3;
NeumannBCInd_u = temp4;