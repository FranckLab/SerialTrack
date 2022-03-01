%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SerialTrack execute main file
% ===================================================
% Dimension:            2D
% Particle rigidity:    hard 
% Tracking mode:        cumulative
% Syn or Exp:           exp  
% Deformation mode:     laser induced cavitation
%
% ===================================================
% Author: Jin Yang, Ph.D.
% Email: jyang526@wisc.edu -or-  aldicdvc@gmail.com 
% Date: 02/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
close all; clear all; clc; clearvars -global
disp('************************************************');
disp('*** Welcome to SerialTrack Particle Tracking ***');
disp('************************************************');
addpath( './function/','./src/');  
 

%% user defined parameters %%%%%

%%%%% Problem dimension and units %%%%%
MPTPara.DIM = 2;   % problem dimension
MPTPara.xstep = 1; % unit: px
MPTPara.tstep = 1; % unit: 1/frame

%%%%% Code mode %%%%%
MPTPara.mode = 'cum'; % {'inc': incremental mode; 
                      %  'cum': cumulative mode; 
                      %  'dbf': double frame}

MPTPara.parType = 'hard'; % {'hard': hard particle; 
                          %  'soft': soft particle}

disp('************************************************');
disp(['Dimention: ',num2str(MPTPara.DIM)]);
disp(['Tracking mode: ',MPTPara.mode]);
disp(['Particle type: ',MPTPara.parType]);
disp('************************************************'); fprintf('\n');

%%%%% Image binary mask file %%%%%
MaskFileLoadingMode = 2; % {0: No mask file
                         %  1: loading only one mask file for all frames; 
                         %  2: loading one mask file for each single frame;
                         %  3: loading matlab mat file for all frames}

if MaskFileLoadingMode == 3
    im_roi_mask_file_path = '.\img_par2track_lung\im_roi.mat';  % TODO: Path for the loaded MATLAB mat file as the mask file
else
    im_roi_mask_file_path = '';  % If there is no mask mat file, leave it as empty;
end

%%%%% Particle detection parameters %%%%%
%%%%% Bead parameters %%%%%
BeadPara.thres = 0.5 ;          % Threshold for detecting particles
BeadPara.beadSize = 3;          % Estimated radius of a single particle
BeadPara.minSize = 2;           % Minimum radius of a single particle
BeadPara.maxSize = 20;          % Maximum radius of a single particle
BeadPara.winSize = [5, 5];      % By default
BeadPara.dccd = [1,1];          % By default
BeadPara.abc = [1,1];           % By default
BeadPara.forloop = 1;           % By default
BeadPara.randNoise = 1e-7;      % By default
BeadPara.PSF = [];              % PSF function; Example: PSF = fspecial('disk', Beadpara.beadSize-1 ); % Disk blur
BeadPara.distMissing = 2;       % Distance threshold to check whether particle has a match or not 
BeadPara.color = 'black';       % Bead color: 'white' -or- 'black'


%% SerialTrack particle tracking

%%%%% Multiple particle tracking (MPT) Parameter %%%%%
MPTPara.f_o_s = 60;              % Size of search field: max(|u|,|v|)
MPTPara.n_neighborsMax = 5;      % Max # of neighboring particles
MPTPara.n_neighborsMin = 1;      % Min # of neighboring particles
MPTPara.locSolver = 1;           % Local solver: 1-topology-based feature; 2-histogram-based feature first and then topology-based feature;
MPTPara.gbSolver = 1;            % Global step solver: 1-moving least square fitting; 2-global regularization; 3-ADMM iterations
MPTPara.smoothness = 1e-3;       % Coefficient of regularization
MPTPara.outlrThres = 1.5;        % Threshold for removing outliers in TPT
MPTPara.maxIterNum = 20;         % Max ADMM iteration number
MPTPara.iterStopThres = 1e-2;    % ADMM iteration stopping threshold
MPTPara.strain_n_neighbors = 20; % # of neighboring particles used in strain gauge
MPTPara.strain_f_o_s = 60;       % Size of virtual strain gauge
MPTPara.usePrevResults = 0;      % Whether use previous results or not: 0-no; 1-yes;  


%%%% Postprocessing: merge trajectory segments %%%%%
distThres = 1;           % distance threshold to connect split trajectory segments
extrapMethod = 'pchip';  % extrapolation scheme to connect split trajectory segments
                         % suggestion: 'nearest' for Brownian motion                          
minTrajSegLength = 10;   % the minimum length of trajectory segment that will be extrapolate 
maxGapTrajSeqLength = 0; % the max frame# gap between connected trajectory segments


%%%%% Execute SerialTrack particle tracking %%%%%
if strcmp(MPTPara.mode,'inc')==1
    run_Serial_MPT_2D_hardpar_inc;
elseif strcmp(MPTPara.mode,'cum')==1
    run_Serial_MPT_2D_hardpar_cum;    
elseif strcmp(MPTPara.mode,'dbf')==1
    run_Serial_MPT_2D_hardpar_dbf;
end
 

%% Other post-processing codes
close all;
addpath('./SLMtools/SLMtools/');

%%%%% Load R vs. t data %%%%%
load('D:\MATLAB\ALTPT2\img_cav_20210318\R_data.mat');
post_yCen = mean(CircleFitPar(end-30:end,1));
post_xCen = mean(CircleFitPar(end-30:end,2));
 

tstep = 1/1e6; % 1/(sampling frequency)  
xstep = 3.2e-6; % m2px ratio 

figure, plot(1e6*[1:1:length(Rnew)]'*tstep, 1e6*Rnew*xstep,'.');
xlabel('Time t(us)'); ylabel('R(um)');
set(gca,'fontsize',18);

R2 = Rnew*xstep;
t = ([1:1:length(Rnew)]-2)*tstep; 
R2= R2(:)'; t2=t(:)';

deltat = tstep/5; % 1e-7 s
t = 0 : deltat : (256-1)*tstep;  % (256-round(RmaxTimeLoc/tstep))*tstep;
R2fit_pp = interp1(t2(1:end), R2(1:end) ,'pchip','pp');

R2_fit = ppval(R2fit_pp,t2);
figure, plot(t2*1e6, R2_fit*1e6,'k','linewidth',1);
xlabel('Time t(us)'); ylabel('R(um)'); set(gca,'fontsize',18);


%%
% RGB color for lines
colorRGBList = [0, .447, .741;  .85, .325, .098;  .929, .694, .125; 
                .494, .184, .556;  .466, .674, .188]; 

tempPlotNum = 0;
            
            
close all; figure,
for ImgSeqNum = [ 4, 17, 32 ]-2
    
    tempPlotNum = tempPlotNum+1;
    
    post_yCen = xstep * mean(CircleFitPar(end-30:end,1));  % Bubble centroid x-coord
    post_xCen = xstep * mean(CircleFitPar(end-30:end,2));  % Bubble centroid y-coord
    % post_yCen = xstep *CircleFitPar(2+ImgSeqNum,1); 
    % post_xCen = xstep *CircleFitPar(2+ImgSeqNum,2);

    post_parCoord_def = xstep*parCoord_prev{ImgSeqNum};
    post_uv_B2A_temp = xstep*uv_B2A_prev{ImgSeqNum-1};
    
    post_u_x = -post_uv_B2A_temp(:,1); post_u_y = -post_uv_B2A_temp(:,2);
    post_x = (post_parCoord_def(:,1)-post_xCen); post_y = (post_parCoord_def(:,2)-post_yCen);
    
    post_r = sqrt( post_x.^2 + post_y.^2 );
    
    post_theta = atan2( post_y, post_x );
    post_u_r = post_u_x.*cos(post_theta) + post_u_y.*sin(post_theta);
    post_u_theta = -post_u_x.*sin(post_theta) + post_u_x.*cos(post_theta);
    
    post_v_r = post_u_r / tstep;
 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1);
    %%%%% velocity %%%%%
    % hold on, plot(1e6*post_r, post_v_r, '.','color',[colorRGBList(tempPlotNum,:)])
    % 
    % [slm,xp,yp] = slmengine(1e6*post_r ,  post_v_r ,'plot','off','concaveup','off','knots', 3 ,'regularization','cross');
    %     hold on; plot(xp,yp,'-','linewidth',2,'color',[colorRGBList(tempPlotNum,:)]);
    
    %%%%% disp %%%%%
    hold on, plot(1e6*post_r, 1e6*post_u_r, '.','color',[colorRGBList(tempPlotNum,:)])
    
    [slm,xp,yp] = slmengine(1e6*post_r,  1e6*post_u_r ,'plot','off','concaveup','off','knots', 3 ,'regularization','cross');
        hold on; plot(xp,yp,'-','linewidth',2,'color',[colorRGBList(tempPlotNum,:)]);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    figure, quiver( 1e6*post_x, 1e6*post_y, 1e6*post_u_x , 1e6*post_u_y ,'color',[colorRGBList(tempPlotNum,:)]);
    hold on; viscircles([0,0]*1e6, Rnew(ImgSeqNum+2)*xstep*1e6,'LineStyle','--','color','k');


end

figure(1);
xlabel('r (um)'); ylabel('ur (um)'); set(gca,'fontsize',18);  axis auto;
axis([40,800,-10,30]); box on;

figure(2);
xlabel('x (um)'); ylabel('y (um)'); set(gca,'fontsize',18);   
axis equal; axis([-500,800,-400,400]); box on; 

figure(3);
xlabel('x (um)'); ylabel('y (um)'); set(gca,'fontsize',18);   
axis equal; axis([-600,800,-400,400]); box on; 

figure(4);
xlabel('x (um)'); ylabel('y (um)'); set(gca,'fontsize',18);   
axis equal; axis([-500,800,-400,400]); box on; 


% figure, quiver( post_x, post_y, post_r, post_r);
% figure, quiver( post_x, post_y, post_theta, post_theta);




%% %%%%% Incremental tracking ratio %%%%%
disp('%%%%% Calculate incremental tracking ratio %%%%%'); fprintf('\n');
track_ratio = zeros(length(Img)-1,1);
DefType = 'exp'; defList = [2:1:length(Img)]';
  
for ImgSeqNum = 2 : length(Img)
    track_A2B = track_A2B_prev{ImgSeqNum-1}; 
    track_ratio(ImgSeqNum-1) = length(track_A2B(track_A2B>0))/size(parCoord_prev{ImgSeqNum},1);      
end
 
fig=figure; ax=axes; hold on; plot(defList,track_ratio,'r^-.','linewidth',1);
adjust_fig(fig,ax,'','',''); box on; title('');
xlabel('Frame #'); ylabel('Incremental tracking ratio');
axis([2,length(file_name),0,1]);


%%%%% Save results %%%%%
disp('%%%%%% ALTPT hard particle tracking: Done! %%%%%%'); fprintf('\n');
results_file_name = 'results_hardpar.mat';
save(results_file_name,'parCoord_prev','uv_B2A_prev','track_A2B_prev','track_B2A_prev');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%% Visualize tracked incremental displacement of each frame %%%%%
disp('%%%%% Plot tracked incremental deformations %%%%%'); fprintf('\n');

%%%%% Experimental parameters %%%%%
try xstep = MPTPara.xstep; catch, xstep = 1; end % unit: um/px
try tstep = MPTPara.tstep; catch, tstep = 1; end % unit: us  
% ImgSeqNum  % Frame #
 
%%%%% Plot tracked incremental displacement field %%%%%
%%%%% Make a video %%%%%
v = VideoWriter('video_2D_inc.avi'); v.FrameRate = 5; open(v); figure,

for ImgSeqNum = 2:length(Img)
    
    % Displacement from tracked particles on deformed frame
    disp_A2B_parCoordB = -uv_B2A_prev{ImgSeqNum-1};
    parCoordB = parCoord_prev{ImgSeqNum};

    %%%%% Plot displacements %%%%%
    clf, plotCone2(parCoordB(:,1)*xstep,parCoordB(:,2)*xstep,disp_A2B_parCoordB(:,1)*xstep/tstep ,disp_A2B_parCoordB(:,2)*xstep/tstep );
    set(gca,'fontsize',18); view(2); box on; axis equal; axis tight; set(gca,'YDir','reverse');
    title(['Tracked velocity (#',num2str(ImgSeqNum),')'],'fontweight','normal');
    xlabel('x'); ylabel('y');
    axis(xstep*[MPTPara.gridxyROIRange.gridx(1), MPTPara.gridxyROIRange.gridx(2), ...
          MPTPara.gridxyROIRange.gridy(1), MPTPara.gridxyROIRange.gridy(2) ]);
     
    frame = getframe(gcf);
    writeVideo(v,frame);
    
end
close(v);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%% Compute trajectory %%%%%
disp('%%%%% Merge tracked trajectory segments %%%%%'); fprintf('\n');

%%%%% Initialization %%%%%
try xstep = MPTPara.xstep; catch xstep = 1; end % um2px ratio
try tstep = MPTPara.tstep; catch tstep = 1; end % time gap between consecutive frames
trajInd = 0; % index of trajectory segments
% distThres = 3; % distance threshold to connect split trajectory segments
% extrapMethod = 'pchip';  % extrapolation scheme to connect split trajectory segments
%                          % suggestion: 'nearest' for Brownian motion                          
% minTrajSegLength = 4;    % the minimum length of trajectory segment that will be extrapolate 
% maxGapTrajSeqLength = 0; % the max frame# gap between connected trajectory segments


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Compute and collect all trajectory segments %%%%%
for tempk = 1 : size(parCoord_prev,1)  % Find trajectories passing through particles in frame #tempk
    
    parCoordCurr = parCoord_prev{tempk}; % Find trajectories passing parCoordCurr
    clear parCoordTrajCurr; parCoordTrajCurr = cell(size(parCoordCurr,1),1); % Initialize a cell structure to store #tempk trajectories
    
    % Add particles in frame #tempk to "parCoordTrajCurr"
    for tempj = 1:length(parCoordCurr)
        parCoordTrajCurr{tempj}(tempk,1:2) = parCoordCurr(tempj,1:2);
    end
    
    for ImgSeqNum = (1+tempk) : size(parCoord_prev,1) % For each later tracked incremental deformation
        
        parCoordB = parCoord_prev{ImgSeqNum}; % Current particle coordinates in frame #ImgSeqNum
        parCoordB_Ind = []; parCoordCurr_Ind = []; % Initialize index
        track_B2A_Curr = track_B2A_prev{ImgSeqNum-1}; % Current tracking in frame #ImgSeqNum
        for tempi = 1:length(track_B2A_Curr)
            try
                % Back propogation
                if ImgSeqNum > 2
                    for tempj = (ImgSeqNum-2) : -1 : tempk
                        track_B2A_Last = track_B2A_prev{tempj};
                        track_B2A_Curr(tempi) = track_B2A_Last(track_B2A_Curr(tempi));
                    end
                end
                % Check whether the trajectory has already been plotted previously
                if tempk>1
                    track_B2A_SecLast = track_B2A_prev{tempk-1};
                    track_B2A_SecLast_temp = track_B2A_SecLast(track_B2A_Curr(tempi)); % if this is not 0, means it's already been plotted
                else % Trajectories from first two frames (tempk=1) will be added
                    track_B2A_SecLast_temp = 0;
                end
                % Assign index values
                if (track_B2A_Curr(tempi)>0) && (track_B2A_SecLast_temp==0)
                    parCoordB_Ind = [parCoordB_Ind; tempi];
                    parCoordCurr_Ind = [parCoordCurr_Ind; track_B2A_Curr(tempi)];
                end
            catch
            end
        end
        
        for tempj = 1:length(parCoordCurr_Ind) % Add found trajectories to cell structure "parCoordTraj"
            parCoordTrajCurr{parCoordCurr_Ind(tempj)}(ImgSeqNum,1:2) = parCoordB(parCoordB_Ind(tempj),1:2);
        end
    end
    
    for parInd = 1:size(parCoordTrajCurr,1)
        wayPoints = parCoordTrajCurr{parInd};
        if ~isempty(wayPoints)
            wayPoints(wayPoints(:,1)==0,:) = wayPoints(wayPoints(:,1)==0,:)*nan;
            wayPoints = [wayPoints; nan(size(parCoord_prev,1)-size(wayPoints,1),2)];
            trajInd = trajInd + 1;
            parCoordTraj{trajInd} = wayPoints;
        end
    end
      
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Merge trajectory segments %%%%%
% Find the starting point and length of each trajectory segment
parCoordTrajPara = []; parCoordTraj = parCoordTraj(:);
for tempi = 1 : size(parCoordTraj,1)
    parCoordTrajCurr = parCoordTraj{tempi};
    parCoordTrajPara(tempi,1) = find(isnan(parCoordTrajCurr(:,1))==0, 1, 'first');
    parCoordTrajPara(tempi,2) = sum(1-isnan(parCoordTrajCurr(:,1)));
end

%%%%% Try to merge trajectory segments by extrapolating the particle position %%%%%
hbar = waitbar(0,'wait');
for tempMergeTime = 1:4 % Try to merge four times
    
    for tempm = 0:maxGapTrajSeqLength % tempm is the # of missing particles between trajectory segments
        
        for tempk = (size(parCoord_prev,1)-1) : -1 : minTrajSegLength % For trajectory segments with length [(size(parCoord_prev,1)-1) : -1 : minTrajSegLength]
            
            [row,col] = find( parCoordTrajPara(:,2)==tempk ); % Find trajectory segments whose "length==tempk"
            [row1,~] = find( parCoordTrajPara(:,2)<size(parCoord_prev,1)+1-tempk & parCoordTrajPara(:,2)>0 ); % Find trajectory setments with length requirement
             
            for tempi = 1:length(row) % For each trajectory segment whose "length==tempk"
                 
                tempWatibar = tempi/length(row)/4/(maxGapTrajSeqLength+1)/((size(parCoord_prev,1)-1)-minTrajSegLength+1) + ...
                    ((size(parCoord_prev,1)-1)-tempk+1)/4/(maxGapTrajSeqLength+1)/((size(parCoord_prev,1)-1)-minTrajSegLength+1) + ...
                    (tempm)/4/(maxGapTrajSeqLength+1) + (tempMergeTime-1)/4;
                
                waitbar(tempWatibar);
                 
                parCoordTrajMat = cell2mat( parCoordTraj );
                parCoordTrajCurr = parCoordTraj{row(tempi)}; % For each trajectory segment whose "length==tempk"
                parCoordTrajCurrPred_x = fillmissing(parCoordTrajCurr(:,1),extrapMethod); % fill missing using 'pchip' method
                parCoordTrajCurrPred_y = fillmissing(parCoordTrajCurr(:,2),extrapMethod); % fill missing using 'pchip' method
                
                % figure, plot(parCoordTrajCurrPred_x,parCoordTrajCurrPred_y,'-o');
                
                
                %%%%% Find all probable trajectory segments in the positive direction %%%%%
                if (sum(parCoordTrajPara(row(tempi),1:2))+tempm<(size(parCoord_prev,1)+1)) && (sum(parCoordTrajPara(row(tempi),1:2))+tempm>0)
                    
                    [row2,col2] = find( parCoordTrajPara(:,1) == sum(parCoordTrajPara(row(tempi),1:2))+tempm ); % starting point requirement
                    row3 = intersect(row1,row2);
                    tempxx = size(parCoord_prev,1)*(row3-1) +(sum(parCoordTrajPara(row(tempi),1:2)))+tempm; tempxx=tempxx'; tempxx=tempxx(:); % Index in parCoordTrajMat
                    tempyy = parCoordTrajMat(tempxx, 1:2);
                    
                    % hold on; plot(tempyy(:,1),tempyy(:,2),'.'); pause;
                    
                    tempzz = sqrt((tempyy(:,1)-parCoordTrajCurrPred_x( sum(parCoordTrajPara(row(tempi),1:2))+tempm )).^2 + ...
                        (tempyy(:,2)-parCoordTrajCurrPred_y( sum(parCoordTrajPara(row(tempi),1:2))+tempm )).^2);
                    
                    [tempzzmin,tempzzminind] = min(tempzz);
                    if tempzzmin < distThres % Find the continuing trajectory segment %JY!!!! threshold distance 3
                        
                        % Merge trajectory segment
                        parCoordTraj{row(tempi)}( parCoordTrajPara(row3(tempzzminind),1) : parCoordTrajPara(row3(tempzzminind),1)+parCoordTrajPara(row3(tempzzminind),2)-1, 1:2 ) = ...
                            parCoordTraj{row3(tempzzminind)}( parCoordTrajPara(row3(tempzzminind),1) : parCoordTrajPara(row3(tempzzminind),1)+parCoordTrajPara(row3(tempzzminind),2)-1, 1:2 );
                        % Remove repeated trajectory segment
                        parCoordTraj{row3(tempzzminind)} = parCoordTraj{row3(tempzzminind)}*nan;
                        
                        % Update varaible "parCoordTrajPara" for trajectory segment length
                        parCoordTrajPara(row(tempi),2) = parCoordTrajPara(row(tempi),2) + parCoordTrajPara(row3(tempzzminind),2) + tempm;
                        parCoordTrajPara(row3(tempzzminind),1:2) = [0,0];
                        
                        %Fillmissing parCoordTraj{row(tempi)}
                        temp = parCoordTraj{row(tempi)}(parCoordTrajPara(row(tempi),1) : sum(parCoordTrajPara(row(tempi),1:2))-1, 1:2);
                        temp_x = fillmissing(temp(:,1),extrapMethod); % fill missing
                        temp_y = fillmissing(temp(:,2),extrapMethod); % fill missing
                        parCoordTraj{row(tempi)}(parCoordTrajPara(row(tempi),1) : sum(parCoordTrajPara(row(tempi),1:2))-1, 1:2) = [temp_x, temp_y];
                        
                    end
                    
                end
                
                
                %%%%% Find all probable trajectory segments in the negative direction %%%%%
                if (sum(parCoordTrajPara(row(tempi),1:2))+tempm<(size(parCoord_prev,1)+1)) && ((parCoordTrajPara(row(tempi),1)-1)-tempm>0)
                    
                    [row2,col2] = find( sum(parCoordTrajPara,2) == parCoordTrajPara(row(tempi),1)-tempm ); % ending point requirement
                    row3 = intersect(row1,row2);
                    tempxx = size(parCoord_prev,1)*(row3-1) + (parCoordTrajPara(row(tempi),1)-1)-tempm; tempxx=tempxx'; tempxx=tempxx(:); % Index in parCoordTrajMat
                    tempyy = parCoordTrajMat(tempxx, 1:2);
                    
                    % hold on; plot(tempyy(:,1),tempyy(:,2),'.');
                    tempzz = sqrt( ( tempyy(:,1)-parCoordTrajCurrPred_x( (parCoordTrajPara(row(tempi),1)-1)-tempm ) ).^2 + ...
                        ( tempyy(:,2)-parCoordTrajCurrPred_y( (parCoordTrajPara(row(tempi),1)-1)-tempm ) ).^2 );
                    
                    [tempzzmin,tempzzminind] = min(tempzz);
                    if tempzzmin < distThres % Find the continuing trajectory segment %JY!!!! threshold distance 3
                        % Merge trajectory segment
                        parCoordTraj{row(tempi)}( parCoordTrajPara(row3(tempzzminind),1) : parCoordTrajPara(row3(tempzzminind),1)+parCoordTrajPara(row3(tempzzminind),2)-1, 1:2 ) = ...
                            parCoordTraj{row3(tempzzminind)}( parCoordTrajPara(row3(tempzzminind),1) : parCoordTrajPara(row3(tempzzminind),1)+parCoordTrajPara(row3(tempzzminind),2)-1, 1:2 );
                        % Remove repeated trajectory segment
                        parCoordTraj{row3(tempzzminind)} = parCoordTraj{row3(tempzzminind)}*nan;
                        % Update varaible "parCoordTrajPara" for both trajectory segment starting point and its length
                        parCoordTrajPara(row(tempi),2) = parCoordTrajPara(row(tempi),2) + parCoordTrajPara(row3(tempzzminind),2) + tempm;
                        parCoordTrajPara(row(tempi),1) = parCoordTrajPara(row3(tempzzminind),1);
                        parCoordTrajPara(row3(tempzzminind),1:2) = [0,0];
                        
                        %Fillmissing parCoordTraj{row(tempi)}
                        temp = parCoordTraj{row(tempi)}(parCoordTrajPara(row(tempi),1) : sum(parCoordTrajPara(row(tempi),1:2))-1, 1:2);
                        temp_x = fillmissing(temp(:,1),extrapMethod); % fill missing
                        temp_y = fillmissing(temp(:,2),extrapMethod); % fill missing
                        parCoordTraj{row(tempi)}(parCoordTrajPara(row(tempi),1) : sum(parCoordTrajPara(row(tempi),1:2))-1, 1:2) = [temp_x, temp_y];
                        
                    end 
                    
                end
                
            end % End of "for each trajectory segment whose "length==tempk" "
            
            
            
        end % End of "for trajectory segments with length [(size(parCoord_prev,1)-1) : -1 : minTrajSegLength]"
        
    end % End of tempm
    
end % End of "for tempMergeTime = 1:5 % Try to merge four times"
close(hbar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plot tracked trajectories %%%%%
disp('%%%%% Plot tracked trajectories %%%%%'); fprintf('\n');

figure,
% CN = [36 63 146; 20 76 156; 28 100 175; 9 115 186; 23 128 196; 33 148 210;
%     25 158 218; 19 172 226; 25 186 229; 69 196 221; 118 205 214; 157 216 215;
%     169 220 217; 193 229 224; 216 237 233; 234 246 245]/255;
% CN = CN(end-size(parCoord_prev,1):end, 1:3);


for tempi = 1:size(parCoordTrajPara,1)
    
    wayPoints = parCoordTraj{tempi};
    
    if (length(wayPoints(isnan(wayPoints(:,1))<1,1))+1)<4
        hold on; line(xstep*wayPoints(isnan(wayPoints(:,1))<1,1),xstep*wayPoints(isnan(wayPoints(:,1))<1,2),'linewidth',1.2); view(2); % straight lines
    else
        hold on; fnplt(cscvn(xstep*wayPoints(isnan(wayPoints(:,1))<1,:)'),'',1.2);
    end
     
    %%%%% Codes to plot trajectories with frame-dependent color %%%%% 
    % if sum(1-isnan(wayPoints(:,1)))>1  % Don't show if there is only one point on the trajectory
    %     hold on; plot(xstep*wayPoints(:,1),xstep*wayPoints(:,2),'r.','markersize',5);
    % end
    % 
    % for tempj = 1:size(parCoord_prev,1)-1
    %     hold on; line(xstep*[wayPoints(tempj,1),wayPoints(tempj+1,1)], ...
    %         xstep*[wayPoints(tempj,2),wayPoints(tempj+1,2)],'linewidth',1.2, 'color', CN(tempj,:) );
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
end
 
set(gca,'fontsize',18); view(2); box on; axis equal; axis tight;  
title('Tracked particle trajectory','fontweight','normal');
set(gca,'YDir','reverse'); xlabel('x'); ylabel('y');
axis(xstep*[MPTPara.gridxyROIRange.gridx(1), MPTPara.gridxyROIRange.gridx(2), ...
      MPTPara.gridxyROIRange.gridy(1), MPTPara.gridxyROIRange.gridy(2) ]);
   

  
%%
%%%%% Compute cumulative tracking ratio from merged trajectories %%%%%
disp('%%%%% Plot tracked cumulative displacements %%%%%'); fprintf('\n');
parCoordTrajMat = cell2mat( parCoordTraj );

[row1,col1] = find(isnan(parCoordTrajMat(1:length(Img):end,1))==0);
trackParCum_ind = row1;
trackParCum_track_ratio = [];

for ImgSeqNum = 2:length(Img)
    [row2,col2] = find(isnan(parCoordTrajMat(ImgSeqNum:length(Img):end,1))==0);
    trackParCum_ind = intersect(row2,trackParCum_ind);
    trackParCum_track_ratio(ImgSeqNum-1) = length(trackParCum_ind) /size(parCoord_prev{1},1);
end

defList = [2:1:length(Img)];
fig=figure; ax=axes; hold on; plot(defList,trackParCum_track_ratio,'bs--','linewidth',1);
adjust_fig(fig,ax,'','',''); box on; title('');
xlabel('Image #'); ylabel('Cumulative tracking ratio');
axis([2,length(Img),0,1]);


%%%%% Plot tracked cumulative displacement field %%%%%
%%%%% Make a video %%%%%
v = VideoWriter('video_2D_inc_cum.avi'); v.FrameRate = 5; open(v); figure,
for ImgSeqNum = 2:length(Img)
    
    parCoordA = parCoordTrajMat(1:length(Img):end,1:2);
    parCoordB = parCoordTrajMat(ImgSeqNum:length(Img):end,1:2);
    parCoordACum = parCoordA(trackParCum_ind,:);
    parCoordBCum = parCoordB(trackParCum_ind,:);
    disp_A2BCum = parCoordBCum - parCoordACum;
    
    % ----- Cone plot grid data: displecement -----
    clf; plotCone2(xstep*parCoordBCum(:,1),xstep*parCoordBCum(:,2),xstep*disp_A2BCum(:,1),xstep*disp_A2BCum(:,2));
    set(gca,'fontsize',18); view(2); box on; axis equal; axis tight; set(gca,'YDir','reverse');
    title(['Tracked cumulative disp (#',num2str(ImgSeqNum),')'],'fontweight','normal');
    xlabel('x'); ylabel('y');
    axis(xstep*[MPTPara.gridxyROIRange.gridx(1), MPTPara.gridxyROIRange.gridx(2), ...
          MPTPara.gridxyROIRange.gridy(1), MPTPara.gridxyROIRange.gridy(2) ]);
     
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Modify codes below to plot interpolated displacements and strains on a uniform grid mesh');
pause; 


ImgSeqNum = 2; % Frame #


%%%%% Previously tracked displacement field %%%%%
parCoordA = parCoordTrajMat(1:length(Img):end,1:2);
parCoordB = parCoordTrajMat(ImgSeqNum:length(Img):end,1:2);
parCoordACum = parCoordA(trackParCum_ind,1:2);
parCoordBCum = parCoordB(trackParCum_ind,1:2);
disp_A2BCum = parCoordBCum - parCoordACum;
  
%%%%% Interpolate scatterred data to gridded data %%%%%
sxy = min([round(0.5*MPTPara.f_o_s),20])*[1,1]; % Step size for griddata
smoothness = 1e-3; % Smoothness for regularization; "smoothness=0" means no regularization

[x_Grid_refB,y_Grid_refB,u_Grid_refB]=funScatter2Grid2D(parCoordBCum(:,1),parCoordBCum(:,2),disp_A2BCum(:,1),sxy,smoothness);
[~,~,v_Grid_refB]=funScatter2Grid2D(parCoordBCum(:,1),parCoordBCum(:,2),disp_A2BCum(:,2),sxy,smoothness);

% Apply ROI image mask
[u_Grid_refB, v_Grid_refB] = funRmROIOutside(x_Grid_refB,y_Grid_refB,MPTPara.ImgRefMask,u_Grid_refB,v_Grid_refB);
% Build a displacement vector
uv_Grid_refB_Vector=[u_Grid_refB(:),v_Grid_refB(:)]'; uv_Grid_refB_Vector=uv_Grid_refB_Vector(:);
% Calculate deformation gradient
D_Grid = funDerivativeOp(size(x_Grid_refB,1),size(x_Grid_refB,2),mean(sxy)); % Central finite difference operator
F_Grid_refB_Vector=D_Grid*uv_Grid_refB_Vector; % {F}={D}{U}


%%%%% Cone plot grid data: displecement %%%%%
figure, plotCone2(xstep*x_Grid_refB,xstep*y_Grid_refB,u_Grid_refB*xstep ,v_Grid_refB*xstep );
set(gca,'fontsize',18); view(2); box on; axis equal; axis tight; set(gca,'YDir','reverse');
title('Tracked cumulative displacement','fontweight','normal');
axis(xstep*[MPTPara.gridxyROIRange.gridx(1), MPTPara.gridxyROIRange.gridx(2), ...
      MPTPara.gridxyROIRange.gridy(1), MPTPara.gridxyROIRange.gridy(2) ]);

  
%%%%% Generate a FE-mesh %%%%%
[coordinatesFEM_refB,elementsFEM_refB] = funMeshSetUp(x_Grid_refB*xstep,y_Grid_refB*xstep);

%%%%% Cone plot grid data: displacement %%%%%
Plotdisp_show(uv_Grid_refB_Vector*xstep , coordinatesFEM_refB, elementsFEM_refB,[],'NoEdgeColor');
 
%%%%% Cone plot grid data: infinitesimal strain %%%%%
Plotstrain_show(F_Grid_refB_Vector, coordinatesFEM_refB, elementsFEM_refB,[],'NoEdgeColor',xstep,tstep);
 

  