function ret= HAR_frcst( Y,X,insmpl )
%This function computes one-step ahead rolling window HAR forecasts with
%re-estimation with the HAR model:
%Y(t)=w+b_d*X_d(t-1)+b_w*X_w(t-1)+b_m*X_m(t-1)+u(t)
%Inputs: Y: a T-by-1 vector used as the l.h.s. of the HAR model - RK 
%        X: a T-by-1 vector used as the r.h.s. of the HAR model
%        insmpl: a integer indicating the size of the initial estimation
%        window - 504
%Outputs: ret: a (T-insmpl)-by-2 vector. The first column is the true value
%and the second column stores the forecasted value
% 2013 x 2 

%Input Checking
if length(Y)~=length(X)
   error('Size mismatch between Y and X');
end
if insmpl>=length(Y)
   error('In-sample period should be shorter than the length of the data'); 
end

%Extract the length of the dataset
T=length(Y); % 2517
%Get the number of re-estimations of the HAR model
N=T-insmpl; % = 2517 - 504 = 2013
%Initialize a vector to compute results
ret=zeros(N,2);% 2013 x 2
%Set the first column of ret to be true values
ret(:,1)=Y(insmpl+1:T);
% første kolonne i ret = Y[504 + 1:2517] = RK[505:2517] = RK[out of sample values]
% 
for i=1:N % 1,2,...,2013 (for each forecasted day)
    %Get temporary Y and X for each re-estimation of HAR. Think about the
    %index i:insmpl+i-1!
    %rolling insmple values 
    Ytemp=Y(i:insmpl+i-1); % RK[1:504]=insmpl,RK[2:505],RK[3:506],...,RK[2013:2516] (all have length: 504 = insmpl size)
    Xtemp=X(i:insmpl+i-1); % RV[1:504],RV[2:505],...,RV[2013:2516] (all have length: 504, which is insmpl size)
    %Compute the r.h.s. and l.h.s. variables using the HARX.m function
    [h,XF]= HARX( Ytemp,Xtemp);
    % h: insmpl x 4 matrix = 482 x 4 = [Y, X_daily, X_weekly, X_monthly]
    % XF:  a 1-by-3 vector that can be used to forecast Y(t+1) %
    % coefficients 
    
    %OLS regression by hand...
    %HAR 
    if i == 1
         har=[ ones(insmpl-22,1) h(:,2:4)]; % 482 x 4  matrix 
     % = [ones(504-22 = 482,1) [X_daily X_weekly X_monthly]
     % = [ har1 har2 har3 har4] 
     % hvor har1 bare er en kolonne av lengde 482 med elementet 1
     
         b=(har'*har)\har'*h(:,1); % beta 
     % = (transpose(har)*har) % transpose(har)* h(:,1)
         ret(i,2)=[1 XF]*b; % forecasted value 1 ...  2013
    else %                                      *    Y
     % b
     % 0.0000
    % 0.0258
    % 0.6119
   % -0.1246
     
     %Compute the one-step ahead forecast and store it
     ret(i,2)=[1 XF]*b; % forecasted value 1 ...  2013
    end 
end
end

