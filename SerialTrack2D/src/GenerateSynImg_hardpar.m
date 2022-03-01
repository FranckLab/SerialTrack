%% Run file to generate both reference & deformed synthetic images
% Two images are synthesized with random dot patterns
%
% %%%%% DEFORM IMAGE: hard particles %%%%%
% Particles are called "HARD" means that shapes of these particles don't  
% change during sample's deformation.
% ----------------------------------------------
%
% ----------------------------------------------
% References
% 
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 02/2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; % clc

%% Generate a synthetic image
 
sigma = [1 1]; % desired standard deviation of bead
sizeI = [512,512]; % desired output volumetric image size.
seedingDensity  = 0.012; % beads per pixel
nBeads = round(seedingDensity*prod(sizeI)); % number of beads

% generate seed points from uniform distribution and scale to image size
%%%%% uniform random %%%%%
% x0{1} =  rand(nBeads,2);
% x0{1} = [x0{1}(:,1)*(sizeI(1)-1) + 1, x0{1}(:,2)*(sizeI(2)-1) + 1];
% I0 = seedBeadsN(sigma,x0{1},sizeI);

%%%%% Poisson Disc sampling %%%%%
[x0] = poissonDisc(sizeI,5,nBeads,1);
I0 = seedBeadsN(sigma,x0,sizeI);
x0temp = x0; clear x0; x0{1} = x0temp;

figure, imshow(uint8(255*I0)'); 
%fig_name = ['img_',num2str(1000+i),'.tif'];
% print(filename,'-dtiff');
%imwrite(uint8(ImgCurr)',fig_name);


%% %%%%% DEFORM IMAGE: hard particles %%%%%

nDeformations = 19; % number of deformation increments +1;   +2 for translations

%%%%% applied translation %%%%%
Trans(1) = 0; Trans(2) = 0;  

%%%%% applied stretch %%%%%
StretchRatio(1) = 1;  
StretchRatio(2) = 1; % sqrt(1/L(1));
  
%%%%% applied rotation %%%%%
RotAng = 180/180*pi;

%%%%% applied simple shear %%%%%
SimpleShear = 0;

I{1} = I0;
u = cell(1,length(2:nDeformations)); 

for i = 2:nDeformations
    fprintf('Deformation step: %i / %i \n',i-1,nDeformations-1)
    
    % step = (i-2)/(nDeformations-2); % Translation
    step = (i-1)/(nDeformations-1); % Non-translation
    
    xC = (sizeI)/2; % center of deformation

    %%%%% Stretch %%%%%
%     u{i-1}(:,1) = step*(StretchRatio(1) - 1)*(x0{1}(:,1) - xC(1)) + step*Trans(1);
%     u{i-1}(:,2) = step*(StretchRatio(2) - 1)*(x0{1}(:,2) - xC(2)) + step*Trans(2);

    %%%%% Rotation %%%%%
   u{i-1}(:,1) = (cos(step*RotAng)-1)*(x0{1}(:,1) - xC(1)) - sin(step*RotAng)*(x0{1}(:,2) - xC(2))  + step*Trans(1);
   u{i-1}(:,2) = sin(step*RotAng)*(x0{1}(:,1) - xC(1)) + (cos(step*RotAng)-1)*(x0{1}(:,2) - xC(2))  + step*Trans(2);
      
    %%%%% Affine %%%%%
%     F11 = 0;  
%     F12 = step*SimpleShear;
%     F21 = 0; 
%     F22 = 0;
%     u{i-1}(:,1) = F11*(x0{1}(:,1) - xC(1)) + F12*(x0{1}(:,2) - xC(2)) + step*Trans(1);
%     u{i-1}(:,2) = F21*(x0{1}(:,1) - xC(1)) + F22*(x0{1}(:,2) - xC(2)) + step*Trans(2);

      
    x1{i-1} = x0{1} + u{i-1};
    I{i} = seedBeadsN(sigma,x1{i-1},sizeI) ;
    
end

for i = 1:length(I)
ImgCurr =  uint8(255*( I{i} + 0.05*randn(sizeI) )) ;
figure, imshow(uint8(ImgCurr)'); 
fig_name = ['img_',num2str(1000+i),'.tif'];
% print(filename,'-dtiff');
imwrite(uint8(ImgCurr)',fig_name);
% save(['vol_rot150_sp_',num2str(1000+i),'.mat'],'vol');
end

save(['imposed_disp','.mat'],'u','x0','x1','nBeads','Trans','StretchRatio','RotAng','SimpleShear','sigma','seedingDensity');
disp('------ Done! ------');

 






