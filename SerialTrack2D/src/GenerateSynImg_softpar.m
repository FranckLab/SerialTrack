%% Run file to generate both reference & deformed synthetic images
% Two images are synthesized with random dot patterns
%
% %%%%% DEFORM IMAGE: soft particles %%%%%
% Particles are called "SOFT" means that shapes of these particles also  
% change with sample's deformation.
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

DefType = 'translation' ;   %{ 'translation' ,  'stretch',   'simpleshear',  'rotation' };
seedingDensity = 0.012;  %{ 0.003, 0.006, 0.012 };

%% Generate a synthetic image

sigma = [1 1]; % desired standard deviation of bead
if strcmp(DefType,'rotation')~=1 % if deformation is not rigid body rotation
    sizeI = [ 512 ,512 ];
else
    sizeI = 2*[ 512 ,512 ]; % desired output volumetric image size.
end
% seedingDensity = 0.003; % beads per pixel
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

ImgCurr = uint8(255*( I0 + 0.05*randn(sizeI) )) ;
figure, imshow(uint8(ImgCurr)'); 
imwrite(ImgCurr','img_1001.tif');


%% %%%%% DEFORM IMAGE: soft particles %%%%%

close all;
file_name = 'img_1001.tif'; % reference image file name
f = imread(file_name); f = double(f); % load referemace image
try mex -O ba_interp2.cpp; catch; end
M = sizeI(1); N = sizeI(2);

g = zeros(M,N); % initialize deformed image
ImgCenX = (sizeI(1))/2; ImgCenY = (sizeI(2))/2; % center of deformation
[xGrid,yGrid] = ndgrid(1:M,1:N);

% number of deformation increments
if strcmp(DefType,'translation')==1
    nDeformations = 41 ; 
elseif strcmp(DefType,'stretch')==1
    nDeformations = 21 ;  
elseif strcmp(DefType,'simpleshear')==1
    nDeformations = 21 ; 
elseif strcmp(DefType,'rotation')==1
    nDeformations = 19 ; 
end

%%%%% applied translation %%%%%
Trans(1) = 0; Trans(2) = 0;  
if strcmp(DefType,'translation')==1, Trans(1)=4; end

%%%%% applied stretch %%%%%
StretchRatio(1) = 1;  
StretchRatio(2) = 1; % sqrt(1/L(1));
if strcmp(DefType,'stretch')==1, StretchRatio(1) = 3; end

%%%%% applied rotation %%%%%
RotAng = 0/180*pi;
if strcmp(DefType,'rotation')==1, RotAng = 180/180*pi; end

%%%%% applied simple shear %%%%%
SimpleShear = 0;
if strcmp(DefType,'simpleshear')==1, SimpleShear = 1; end

 
for i = 2:nDeformations
    
    fprintf('Deformation step: %i / %i \n',i-1,nDeformations-1)
    
    step = (i-1)/(nDeformations-1);
    
    if strcmp(DefType,'translation')==1 %%%%% Translation %%%%% 
        F11Inv = 0; F12Inv = 0; F21Inv = 0; F22Inv = 0;  
    end
    if strcmp(DefType,'stretch')==1 %%%%% Stretch %%%%%
        F11Inv = 1/(1+step*(StretchRatio(1)-1))-1; F12Inv = 0; F21Inv = 0; F22Inv = 0;
    end
    if strcmp(DefType,'simpleshear')==1 %%%%% Simple shear %%%%%
    	F11Inv = 0; F12Inv = -step*SimpleShear; F21Inv = 0; F22Inv = 0;
    end
    if strcmp(DefType,'rotation')==1 %%%%% Rotation %%%%%
        F11Inv = (cos(step*RotAng)-1); F12Inv = sin(step*RotAng);
        F21Inv = -sin(step*RotAng); F22Inv = (cos(step*RotAng)-1);
    end
    
    %%%%% Compute deformed x- and y-coordinates %%%%%
    tempxGrid = (1+F11Inv)*(xGrid-ImgCenX) + (F12Inv)*(yGrid-ImgCenX) + ImgCenX + step*(-Trans(1));
    tempyGrid = (F21Inv)*(xGrid-ImgCenX) + (1+F22Inv)*(yGrid-ImgCenX) + ImgCenY + step*(-Trans(2));
     
    %%%%% Store other information %%%%%
    x0{i-1} = [xGrid(:),yGrid(:)];
    x1{i-1} = [tempxGrid(:),tempyGrid(:)];
    u{i-1} = -x1{i-1} + x0{i-1};
    
    g = ba_interp2(f,tempxGrid,tempyGrid,'cubic') + 255*0.05*randn(sizeI);
    if strcmp(DefType,'rotation')==1 %%%%% if Rotation %%%%%
        g = g( round(0.25*M)+1:round(0.75*M), round(0.25*N)+1:round(0.75*N) );
        xGrid1 = xGrid( round(0.25*M)+1:round(0.75*M), round(0.25*N)+1:round(0.75*N) );
        yGrid1 = yGrid( round(0.25*M)+1:round(0.75*M), round(0.25*N)+1:round(0.75*N) ) ;
        xGrid1_000 = xGrid1(1,1); yGrid1_000 = yGrid1(1,1);
        xGrid1 = xGrid1 - xGrid1_000;
        yGrid1 = yGrid1 - yGrid1_000;
        tempxGrid1 = tempxGrid( round(0.25*M)+1:round(0.75*M), round(0.25*N)+1:round(0.75*N) );
        tempxGrid1 = tempxGrid1 - xGrid1_000;
        tempyGrid1 = tempyGrid( round(0.25*M)+1:round(0.75*M), round(0.25*N)+1:round(0.75*N) );
        tempyGrid1 = tempyGrid1 - yGrid1_000;
        x0{i-1} = [xGrid1(:),yGrid1(:)];
        x1{i-1} = [tempxGrid1(:),tempyGrid1(:)];
        u{i-1} = -x1{i-1} + x0{i-1}; 
        
        % pause;
        
    end
    figure, imshow(uint8(g)');
    % figure, plotCone2( x0{1}(:,1), x0{1}(:,2), u{1}(:,1), u{1}(:,2) );
    %%%%% Save figures %%%%%
    if i<10, fig_name = ['img_100',num2str(i),'.tif'];
    elseif i<100,    fig_name = ['img_10',num2str(i),'.tif'];
    elseif i<1000,   fig_name = ['img_1',num2str(i),'.tif'];
    end
    imwrite(uint8(g)',fig_name);
 
     
    
end

save(['imposed_disp','.mat'],'u','x0','x1','nBeads','Trans','StretchRatio','RotAng','SimpleShear','sigma','seedingDensity');
disp('------ Done! ------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(DefType,'rotation')==1  %%%%% if Rotation %%%%%
    close all;
    file_name = 'img_1001.tif'; % reference image file name
    f = imread(file_name); 
    f = f( round(0.25*M)+1:round(0.75*M), round(0.25*N)+1:round(0.75*N) );
    figure, imshow(uint8(f) );
    imwrite(uint8(f) ,'img_1001.tif');
end


