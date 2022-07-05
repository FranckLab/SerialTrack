%% Run file to generate both reference & deformed synthetic volumetric images
% Two volumetic images are synthesized with random dot patterns

% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 02/2020.
% ==============================================

clear; close all; % clc

DefType = 'stretch';
seedingDensity = 300*1e-6; % {10,100,300,1000}*1e-6; 

%% Mohak's code to generate a synthetic volumetic image
 
sigma = [1 1 1]; % desired standard deviation of bead
sizeI = [192,192,192]; % desired output volumetric image size.
% seedingDensity  = 0.000006; % beads per pixel
nBeads = round(seedingDensity*prod(sizeI)); % number of beads

% generate seed points from uniform distribution and scale to image size
%%%%% uniform random %%%%%
% x0{1} =  rand(nBeads,3);
% x0{1} = [x0{1}(:,1)*(sizeI(1)-1) + 1, x0{1}(:,2)*(sizeI(2)-1) + 1,x0{1}(:,3)*(sizeI(3)-1) + 1];
% I0 = seedBeadsN(sigma,x0{1},sizeI);

%%%%% Poisson Disc sampling %%%%%
[x0] = poissonDisc(sizeI,5,nBeads,1);
I0 = seedBeadsN(sigma,x0,sizeI);
x0temp = x0; clear x0; x0{1} = x0temp;
% figure, plot3(x0temp(:,1),x0temp(:,2),x0temp(:,3),'.');
% figure, imagesc3(I0)

%% DEFORM IMAGE

% number of deformation increments
if strcmp(DefType,'translation')==1
    nDeformations = 41 ; 
elseif strcmp(DefType,'stretch')==1
    nDeformations = 41 ;  
elseif strcmp(DefType,'simpleshear')==1
    nDeformations = 21 ; 
elseif strcmp(DefType,'rotation')==1
    nDeformations = 19 ; 
end
 
%%%%% applied translation %%%%%
Trans(1) = 0; Trans(2) = 0; Trans(3) = 0;
if strcmp(DefType,'translation')==1, Trans(1)=4; end

%%%%% applied stretch %%%%%
StretchRatio(1) = 1.0;  
StretchRatio(2) = 1; % sqrt(1/L(1));
StretchRatio(3) = StretchRatio(2);
if strcmp(DefType,'stretch')==1, StretchRatio(1) = 3; end

%%%%% applied rotation %%%%%
RotAng = 0/180*pi;
if strcmp(DefType,'rotation')==1, RotAng = 180/180*pi; end

%%%%% applied simple shear %%%%%
SimpleShear = 0;
if strcmp(DefType,'simpleshear')==1, SimpleShear = 1; end


I{1} = I0;
u = cell(1,length(2:nDeformations)); 

for i = 2:nDeformations
    
    fprintf('Deformation step: %i / %i \n',i-1,nDeformations-1)
    
    % step = (i-2)/(nDeformations-2); % Translation
    step = (i-1)/(nDeformations-1); % Non-translation
    
    xC = (sizeI)/2; % center of deformation

    if strcmp(DefType,'translation')==1 ... %%%%% Translation %%%%%
           ||  strcmp(DefType,'stretch')==1 %%%%% Stretch %%%%%
        u{i-1}(:,1) = step*(StretchRatio(1) - 1)*(x0{1}(:,1) - xC(1)) + step*Trans(1);
        u{i-1}(:,2) = step*(StretchRatio(2) - 1)*(x0{1}(:,2) - xC(2)) + step*Trans(2);
        u{i-1}(:,3) = step*(StretchRatio(3) - 1)*(x0{1}(:,3) - xC(3)) + step*Trans(3);
    elseif strcmp(DefType,'rotation')==1 %%%%% Rotation %%%%%
        u{i-1}(:,1) = (cos(step*RotAng)-1)*(x0{1}(:,1) - xC(1)) - sin(step*RotAng)*(x0{1}(:,2) - xC(2))  + step*Trans(1);
        u{i-1}(:,2) = sin(step*RotAng)*(x0{1}(:,1) - xC(1)) + (cos(step*RotAng)-1)*(x0{1}(:,2) - xC(2))  + step*Trans(2);
        u{i-1}(:,3) = step*(StretchRatio(3) - 1)*(x0{1}(:,3) - xC(3)) + step*Trans(3);
    elseif strcmp(DefType,'simpleshear')==1 %%%%% Simple shear %%%%%
        F11 = 0; F12 = step*SimpleShear; F21 = 0; F22 = 0;
        u{i-1}(:,1) = F11*(x0{1}(:,1) - xC(1)) + F12*(x0{1}(:,2) - xC(2)) + step*Trans(1);
        u{i-1}(:,2) = F21*(x0{1}(:,1) - xC(1)) + F22*(x0{1}(:,2) - xC(2)) + step*Trans(2);
        u{i-1}(:,3) = step*(StretchRatio(3) - 1)*(x0{1}(:,3) - xC(3)) + step*Trans(3);
    end
    
    x1{i-1} = x0{1} + u{i-1};
    I{i} = seedBeadsN(sigma,x1{i-1},sizeI) ;
    
end

for i = 1:length(I)
vol{1} =  uint8(255*( I{i} + 0.05*randn(sizeI) )) ;
save(['vol_',num2str(1000+i),'.mat'],'vol');
end

save(['imposed_disp','.mat'],'u','x0','x1','nBeads','Trans','StretchRatio','sigma','seedingDensity');
disp('------ Done! ------');

% %%%%% Plot %%%%% 
% figure, imagesc3(vol{1});




