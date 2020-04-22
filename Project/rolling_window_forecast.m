function [y, ymse ] = rolling_window_forecast( mdl,insmpl,outofsmpl,option )
% [y1, ymse1 ] = rolling_window_forecast( mdl,insmpl,outofsmpl,0 );

%This function computes rolling window one-step-ahead forecast for a given
%ARMA model.
%Inputs: mdl: an estimated ARMA model structure
%     insmpl: a vector of data used to estimate the ARMA model - price
%     level
%  outofsmpl: a vector of data to be predicted - price level
%Outputs: y: a vector of the forecasts obtained for the out-of-sample
%period
%      ymse: estimated mean squared error of the forecasts
%option: option=1 for re-estimation of the  ARMA parameter, and 0 otherwise

%Size of the in-sample period
T=length(insmpl); % 6783
%Size of the out-of-sample period
T1=length(outofsmpl); % 441
%initializing y and ymse
y=zeros(T1,1); % vector of length 441 to store forecasted y (price level)
ymse=zeros(T1,1); % vector to stor ymse of y 
%reconstruct the data vector
data=[insmpl;outofsmpl]; % 7224 1 som inneholder (som rader insmpl data + outofsmpl data]
% med andre ord inneholder data y (price level) for hele period = length
% insmpl + length outofsmpl

%compute forecasts based on rolling window...
for i=1:T1 % 441
    if option==1 % reestimate ARMA parameters for insmpl data 
        mtemp=estimate(mdl,insmpl,'Display','off'); % to reestimate ARMA 
        [y(i),ymse(i)]=forecast(mtemp,1,'Y0',data(1+i:T+i));
    else
        if i==1 % first forecast 
            mtemp=estimate(mdl,insmpl,'Display','off'); % estimated model for insmpl period
            [y(i),ymse(i)]=forecast(mtemp,1,'Y0',data(1+i:T+i)); 
            % forecasted model (1) = 
            % Forecast() 
            % IN:
            % 1 = nr of periods to be forecasted = 1 day
            % data = data[2:T+1] 
            % out: y(1) = Y_(T + 1)
        else
            [y(i),ymse(i)]=forecast(mtemp,1,'Y0',data(1+i:T+i));
        end
    end    
end

