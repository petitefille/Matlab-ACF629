
function [ret,XF]= HARX( Y,X )
%This function is used to compute the l.h.s. and r.h.s. variables for a HAR
%model: Y(t)=w+b_d*X_d(t-1)+b_w*X_w(t-1)+b_m*X_m(t-1)+u(t)
%Input: Y: a T-by-1 vector for the l.h.s. of HAR model = RK, T = 504 =
                     %insmpl size so Y: 504 x 1 matrix, the insmple rolling values
%       X: a T-by-1 vector for the r.h.s. of HAR model = RV which is also 
         % a 504 x 1 matrix, the insmple rolling values 
%Output: ret: a (T-22)-by-4 matrix. [Y, X_d, X_w, X_m] => a (504-22) x 22 matrix
%        = (482 x 22) matrix 
%        XF:  a 1-by-3 vector that can be used to forecast Y(t+1)

%Input Checking
if length(Y)~=length(X) % must be of same size
   error('Size mismatch between Y and X');
end

%Number of rows in Y
r=length(Y); % length of insample rolling window = 504
%Initialize the return matrix for forecasts
ret=zeros(r-22,4); % size: (504-22) x 4 = 482 x 4
% ret = [Y X_d, X_w X_m]

%Retrieve the vector of Y as l.h.s. of HAR
ret(:,1)=Y(23:r); % first column in ret = forecasts 
% RK[23:504] = y

%Retrieve the lagged vector of X as X_daily
ret(:,2)=X(22:r-1); % 2nd column 
% X_d = RV(22:504-1) = RV(22:503)

%Computer weekly and monthly moving averages
Xw=movmean(X,[4,0]); % 504 x 1 
Xm=movmean(X,[21,0]); % 504 x 1
%Store the corresponding X_w and X_m
ret(:,3)=Xw(22:r-1);% = Xw(22:503) - 3rd column => lengde 482 
ret(:,4)=Xm(22:r-1);% = Xm(22:504) - 4rth column => lengde 482 
%Store X_d(T) X_w(T) and X_m(T) as XF. 
XF=[X(r) Xw(r) Xm(r)]; %  XF = [X(504) Xw(504) Xm(504)]
end

% XF = [X(r) Xw(r) Xm(r)]

% XF =

   % 1.0e-03 *

    % 0.1533    0.0564    0.0717
