function [tempu,tempv,tempw] = funInitGuess3(parCoordBPrev,u_B2A_prev,parCoordBCurr,ImgSeqNum) 

% %%%%% Store as previous results
% parCoordBPrev{ImgSeqNum-1} = parCoordB;
% u_B2A_prev{ImgSeqNum-1} = u_B2A_curr;

if ImgSeqNum==3 
     
    FxPrev = scatteredInterpolant([parCoordBPrev{ImgSeqNum-2}(:,1),parCoordBPrev{ImgSeqNum-2}(:,2),parCoordBPrev{ImgSeqNum-2}(:,3)],u_B2A_prev{ImgSeqNum-2}(:,1),'linear','linear');  
    FyPrev = scatteredInterpolant([parCoordBPrev{ImgSeqNum-2}(:,1),parCoordBPrev{ImgSeqNum-2}(:,2),parCoordBPrev{ImgSeqNum-2}(:,3)],u_B2A_prev{ImgSeqNum-2}(:,2),'linear','linear'); 
    FzPrev = scatteredInterpolant([parCoordBPrev{ImgSeqNum-2}(:,1),parCoordBPrev{ImgSeqNum-2}(:,2),parCoordBPrev{ImgSeqNum-2}(:,3)],u_B2A_prev{ImgSeqNum-2}(:,3),'linear','linear'); 
    tempu = 2*FxPrev(parCoordBCurr);
    tempv = 2*FyPrev(parCoordBCurr);
    tempw = 2*FzPrev(parCoordBCurr);
 
elseif ImgSeqNum==4   || ImgSeqNum==5 || ImgSeqNum==6
    
    FxPrev2 = scatteredInterpolant([parCoordBPrev{ImgSeqNum-2}(:,1),parCoordBPrev{ImgSeqNum-2}(:,2),parCoordBPrev{ImgSeqNum-2}(:,3)],u_B2A_prev{ImgSeqNum-2}(:,1),'linear','linear');  
    FyPrev2 = scatteredInterpolant([parCoordBPrev{ImgSeqNum-2}(:,1),parCoordBPrev{ImgSeqNum-2}(:,2),parCoordBPrev{ImgSeqNum-2}(:,3)],u_B2A_prev{ImgSeqNum-2}(:,2),'linear','linear'); 
    FzPrev2 = scatteredInterpolant([parCoordBPrev{ImgSeqNum-2}(:,1),parCoordBPrev{ImgSeqNum-2}(:,2),parCoordBPrev{ImgSeqNum-2}(:,3)],u_B2A_prev{ImgSeqNum-2}(:,3),'linear','linear'); 
    tempu2 = FxPrev2(parCoordBCurr);
    tempv2 = FyPrev2(parCoordBCurr);
    tempw2 = FzPrev2(parCoordBCurr);
    
    FxPrev3 = scatteredInterpolant([parCoordBPrev{ImgSeqNum-3}(:,1),parCoordBPrev{ImgSeqNum-3}(:,2),parCoordBPrev{ImgSeqNum-3}(:,3)],u_B2A_prev{ImgSeqNum-3}(:,1),'linear','linear');  
    FyPrev3 = scatteredInterpolant([parCoordBPrev{ImgSeqNum-3}(:,1),parCoordBPrev{ImgSeqNum-3}(:,2),parCoordBPrev{ImgSeqNum-3}(:,3)],u_B2A_prev{ImgSeqNum-3}(:,2),'linear','linear'); 
    FzPrev3 = scatteredInterpolant([parCoordBPrev{ImgSeqNum-3}(:,1),parCoordBPrev{ImgSeqNum-3}(:,2),parCoordBPrev{ImgSeqNum-3}(:,3)],u_B2A_prev{ImgSeqNum-3}(:,3),'linear','linear'); 
    tempu3 = FxPrev3(parCoordBCurr);
    tempv3 = FyPrev3(parCoordBCurr);
    tempw3 = FzPrev3(parCoordBCurr);
    
    tempu = 2*tempu2 - tempu3; 
    tempv = 2*tempv2 - tempv3;
    tempw = 2*tempw2 - tempw3; 
    
    
elseif ImgSeqNum>6 % When ImgSeqNum > 6: POD predicts next disp U0 from previous results of (ImgSeqNum+[-5:1:-1])
      
    nTime = 5; np = size(parCoordBCurr,1); % "nTime" value 5 is an empirical value, can be changed.
    T_data_u = zeros(nTime,np); T_data_v = zeros(nTime,np); T_data_w = zeros(nTime,np); 
    for tempi = 1:nTime
        FxPrev = scatteredInterpolant([parCoordBPrev{ImgSeqNum-7+tempi}(:,1),parCoordBPrev{ImgSeqNum-7+tempi}(:,2),parCoordBPrev{ImgSeqNum-7+tempi}(:,3)],u_B2A_prev{ImgSeqNum-7+tempi}(:,1),'linear','linear');  
        FyPrev = scatteredInterpolant([parCoordBPrev{ImgSeqNum-7+tempi}(:,1),parCoordBPrev{ImgSeqNum-7+tempi}(:,2),parCoordBPrev{ImgSeqNum-7+tempi}(:,3)],u_B2A_prev{ImgSeqNum-7+tempi}(:,2),'linear','linear');  
        FzPrev = scatteredInterpolant([parCoordBPrev{ImgSeqNum-7+tempi}(:,1),parCoordBPrev{ImgSeqNum-7+tempi}(:,2),parCoordBPrev{ImgSeqNum-7+tempi}(:,3)],u_B2A_prev{ImgSeqNum-7+tempi}(:,3),'linear','linear');  
        T_data_u(tempi,:) = FxPrev(parCoordBCurr);
        T_data_v(tempi,:) = FyPrev(parCoordBCurr);
        T_data_w(tempi,:) = FzPrev(parCoordBCurr);
    end
    nB = 3; t_train = [ImgSeqNum-1-nTime:ImgSeqNum-2]'; t_pre = [ImgSeqNum-1]';
    [u_pred,~,~,~] = funPOR_GPR(T_data_u,t_train,t_pre,nB);
    [v_pred,~,~,~] = funPOR_GPR(T_data_v,t_train,t_pre,nB);
    [w_pred,~,~,~] = funPOR_GPR(T_data_w,t_train,t_pre,nB);
    tempu = u_pred(1,:)'; 
    tempv = v_pred(1,:)'; 
    tempw = w_pred(1,:)';
     
else
    
end
