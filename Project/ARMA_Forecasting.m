%Set up and load the dataset
%clear the pre-existing variables
clearvars;
%set working directory (Use your own path instead! Or browse to where the
%data is stored)
cd('C:\Users\Eier\Documents\Lancaster\ACF609\Lab1');
%load an excel file into MATLAB as table
dmat=readtable('StockIndexFXDATA.xlsx');
data=(dmat{:,13}); % of type double of size 7224 x 1
T=length(data);
%Set an insample period (in this case we use data from 1990 to 2015):
% data goes from 02.11.1990 - 08.09.2017
insmpl=data(year(dmat.DATE)<=2015); % only contains data for column 13 between 2009 - 2015
%Similarly, set an out-of-sample period
outofsmpl=data(year(dmat.DATE)>2015); % dmat.DATE er kolonnnen med navn DATE som inneholder alle datoene
% dmat.DATE(1) 02.01.1990
% year(dmat.DATE(1)) 1990

cd('C:\Users\Eier\Documents\Lancaster\ACF609\Lab4');
%Choose the best model for the insmple period using the BIC. Note that I
% use a function by modifying the Box_Jenkins_Analysis code.
 [pbest,qbest ]= ARMA_optimal( insmpl, 2, 2, 'BIC' ); % assuming cov stationarity
 % [2,1]
 
%Construct the optimized ARMA model
mdl=arima(pbest,0,qbest); % Create univariate autoregressive integrated moving average (ARIMA) model

[estmdl,EstParamCov,logL,info]=estimate(mdl,insmpl); % Estimate ARIMA or ARIMAX model parameters
% estmd: Model containing parameter estimates
% EstParamCov: Variance-covariance matrix of maximum likelihood estimates
% logL: Optimized loglikelihood objective function value
% info: Summary information

 
%    ARIMA(2,0,1) Model (Gaussian Distribution):
 
%                 Value      StandardError    TStatistic      PValue   
%                ________    _____________    __________    ___________

%    Constant    0.048383      0.010597         4.5656        4.981e-06
%    AR{1}         1.6586      0.017968         92.306                0
%    AR{2}       -0.66102      0.017673        -37.403      3.5105e-306
%    MA{1}       -0.78833      0.015177        -51.941                0
%    Variance      2.2027      0.012823         171.78                0



%Compute forecasts of the out-of-sample period based on 1990-2015 data and
%the estimated model. Note that this forecast is the fixed window forecast.
[y, ymse, v]=forecast(estmdl,length(outofsmpl),'Y0',insmpl); % Forecast ARIMA or ARIMAX model responses or conditional variances
% Input: estmdl (estimated model for insample period)
% Output: length(out of sample period)
% y: Minimum mean square error forecasts of response data => prediction
% ymse:  Mean square errors of forecasted responses
% v: Minimum MSE forecasts of conditional variances of future model innovations

%Plot the forecasts with the actual value and 95% confidence bounds
figure
h1 = plot(outofsmpl,'Color',[.3,.3,.3]); % true data for out of sample period
hold on
h2 = plot(1:length(outofsmpl),y,'b','LineWidth',2); % estimated y for out of sample period
h3 = plot(1:length(outofsmpl),y + 1.96*sqrt(ymse),'r:',...
		'LineWidth',2); % upper 95 % CI     
plot(1:length(outofsmpl),y - 1.96*sqrt(ymse),'r:','LineWidth',2); % lower 95 % CI
legend([h1 h2 h3],'Observed','Forecast',...
		'95% Confidence Interval','Location','NorthWest');
title(['Two-Year Forecasts and Approximate 95% '...
			'Confidence Intervals'])
hold off

%Now....do a rolling window one-step-ahead forecast with two competing ARMA models, and
%compare their forecasting performances. Note that we only estimate the
%parameters once in the in-sample period. 

%We use the previous model as model 1, and set our competing model to be,
%for example, an ARMA(1,1) model.
mdl2=arima(1,0,1);

%Write a function for rolling window forecasting...
%Pay attention to the fourth input of the function. 
 [y1, ymse1 ] = rolling_window_forecast( mdl,insmpl,outofsmpl,0 );
 
 [y2, ymse2 ] = rolling_window_forecast( mdl2,insmpl,outofsmpl,0 );

%Compute the forecasting error for the two models:
fe1=y1-outofsmpl;
fe2=y2-outofsmpl;
%Compare the forecasts by computing (1) RMSE (2) MAE, and conduct a
%modified-DM test on the MSE.
RMSE1=sqrt(mean(fe1.^2)); % RMSEP =0.2270
RMSE2=sqrt(mean(fe2.^2)); % 0.1731
MAE1=mean(abs(fe1)); % MAPE 0.1760
MAE2=mean(abs(fe2)); % 0.1414
disp('RMSE and MAE for model 1'); % RMSE and MAE for model 1
disp([RMSE1 MAE1]); %   0.2270    0.1760
disp('RMSE and MAE for model 2'); % RMSE and MAE for model 2
disp([RMSE2 MAE2]); %   0.1731    0.1414
%Retrieve the DM test statistic and the p-value. How do you interpret the
%test results?
[tstat,pv]=dmtest_modified(fe1,fe2);
% tstat = 6.5561 > 0
% pvalue = 1.5508e-10
