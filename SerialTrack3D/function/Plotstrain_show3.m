function Plotstrain_show3(F,coordinatesFEM,elementsFEM,varargin)
%FUNCTION Plotstrain_show(F,coordinatesFEM,elementsFEM)
% To plot DVC solved 3D strain components
% ----------------------------------------------
%
%   INPUT: F                 Deformation gradient tensor: 
%                            F = [F11_node1, F21_node1, F31_node1, ...
%                                 F12_node1, F22_node1, F32_node1, ...
%                                 F13_node1, F23_node1, F33_node1, ...
%                                 ... , ...
%                                 F11_nodeN, F21_nodeN, F31_nodeN, ...
%                                 F12_nodeN, F22_nodeN, F32_nodeN, ...
%                                 F13_nodeN, F23_nodeN, F33_nodeN]';
%          coordinatesFEM    FE mesh coordinates
%          elementsFEM       FE mesh elements
%          DICpara           chosen DVC parameters
%          EdgeColorOrNot    show edge color or not

%   OUTPUT: Plots of exx, eyy, and exy strain fields.
%
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last date modified: 2020.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Generate model
dvcVOI = createpde(1);

% Apply mesh
DT = delaunayTriangulation(coordinatesFEM(:,1), coordinatesFEM(:,2), coordinatesFEM(:,3)); 
geometryFromMesh(dvcVOI,DT.Points',DT.ConnectivityList');
% ------ FEMesh has structure ------
% FEMesh with properties:
% 
%              Nodes: [3x10003 double]
%           Elements: [10x5774 double]
%     MaxElementSize: 9.7980
%     MinElementSize: 4.8990
%      MeshGradation: 1.5000
%     GeometricOrder: 'quadratic'
% ----------------------------------

% plot dvcZOI domain
% ---------------------------------
figure,
subplot(2,3,1), pdeplot3D(dvcVOI,'ColorMapData',F(1:9:end),'FaceAlpha',0.5);

title('Strain $e_{11}$','fontweight','normal','Interpreter','latex');  
set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
set(gcf,'color','w'); colormap jet; colorbar; box on
 
% ---------------------------------
subplot(2,3,2), pdeplot3D(dvcVOI,'ColorMapData',F(5:9:end),'FaceAlpha',0.5);

title('Strain $e_{22}$','fontweight','normal','Interpreter','latex');  
set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
set(gcf,'color','w'); colormap jet; colorbar; box on

% ---------------------------------
subplot(2,3,3), pdeplot3D(dvcVOI,'ColorMapData',F(9:9:end),'FaceAlpha',0.5);

title('Strain $e_{33}$','fontweight','normal','Interpreter','latex');  
set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
set(gcf,'color','w'); colormap jet; colorbar; box on

% ---------------------------------
subplot(2,3,4), pdeplot3D(dvcVOI,'ColorMapData',0.5*(F(2:9:end)+F(4:9:end)),'FaceAlpha',0.5);

title('Strain $e_{12}$','fontweight','normal','Interpreter','latex');  
set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
set(gcf,'color','w'); colormap jet; colorbar; box on

% ---------------------------------
subplot(2,3,5), pdeplot3D(dvcVOI,'ColorMapData',0.5*(F(3:9:end)+F(7:9:end)),'FaceAlpha',0.5);

title('Strain $e_{13}$','fontweight','normal','Interpreter','latex');  
set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
set(gcf,'color','w'); colormap jet; colorbar; box on

% ---------------------------------
subplot(2,3,6), pdeplot3D(dvcVOI,'ColorMapData',0.5*(F(6:9:end)+F(8:9:end)),'FaceAlpha',0.5);

title('Strain $e_{23}$','fontweight','normal','Interpreter','latex');  
set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
set(gcf,'color','w'); colormap jet; colorbar; box on




