% Yt straight but sigma varying
% c linear increase decrease
% drift polynomial growth
% ACF does decrease slowly over time but still significant (then
% hypothesis test will change between rejecting and not rejecting

%Set up and load the dataset
%clear the pre-existing variables
clearvars;
%set working directory (Use your own path instead! Or browse to where the
%data is stored)
%cd('C:\Users\Eier\Documents\Lancaster\ACF609\Lab4')
%load an excel file into MATLAB as table
dmat=readtable('StockIndexFXDATA.xlsx'); % of type table (class(dmat))
data=(dmat{:,13}); % size(data):  7224 x 1 !!!!!! = vector

%Set the significance level of the test
alpha=[0.1,0.05,0.01];

%Augmented Dickey-Fuller test with automatic lag selection based on
%information criterion.Read the ADF_auto.m function description!
sel_method='BIC';
dfout=ADF_auto(data,0,30,sel_method); % dif = 0
%alternatively, you can use a SIC fixed lag selection: 12*floor((T/100)^0.25),
%where T is the number of observations
% S = size(data) = 7224 x 1
% T = S(1) = 7224
% 12*floor((T/100)^0.25)= 24
%  6×3 table
%                   None       Drift      Trend 
%                 ________    _______    _______

%    Test Stat     -2.0757    -5.3246    -5.3251   test statistics
%    10% C.V.      -1.6173    -2.5672    -3.1279   0.01 critical value
%    5% C.V.       -1.9416    -2.8622     -3.413   0.05 critical value
%    1% C.V.       -2.5672    -3.4329    -3.9623   0.1 critical value
%    P-value      0.036608      0.001      0.001
%    Lag used           10         10         10   = optimal nr of lags
%                                                    chosen

% 
% dfout=ADF_auto(data,1,30,sel_method);
% dfout

% dfout =

%   6×3 table

%                  None       Drift      Trend 
                 _______    _______    _______

%    Test Stat    -28.674    -28.673    -28.671
%    10% C.V.     -1.6173    -2.5672    -3.1279
%    5% C.V.      -1.9416    -2.8622     -3.413
%    1% C.V.      -2.5672    -3.4329    -3.9623
%    P-value        0.001      0.001      0.001
%    Lag used           9          9          9
    
    
    
    
%For PP and KPSS test, the selection of lags is more involved. Here I used
%a fixed bandwith for the Newey-West type variance-covariance estimates,
%which may not be the best selection method for the two type of tests. 
Bandwidth=floor(4* (length(data)^(2/9)) ); % 28 and floor(4.9) = 4
%Philip Perron test
%Perform three sets of tests for three different model specifications
[h1,pValue1,stat1,cValue1,reg1]=pptest(data,'Model','AR', 'alpha', alpha,'lags', Bandwidth);
[h2,pValue2,stat2,cValue2,reg2]=pptest(data,'Model','ARD', 'alpha',alpha, 'lags', Bandwidth);
[h3,pValue3,stat3,cValue3,reg3]=pptest(data,'Model','TS',  'alpha', alpha,'lags', Bandwidth);
%Prepare a table to present the test outputs
pptable=[stat1(1) cValue1 pValue1(1); stat2(1) cValue2 pValue2(1) ;stat3(1) cValue3 pValue3(1)];
ppout=array2table(pptable','RowNames', {'Test Stat' ' 10% C.V.' '5% C.V.' '1% C.V.' 'P-value' },...
'VariableNames', {'None' 'Drift' 'Trend'}); 

 %                 None       Drift      Trend 
                 _______    _______    _______

  %  Test Stat    -2.1144    -6.6476    -6.6493
  %  10% C.V.     -1.6173    -2.5672    -3.1279
  %  5% C.V.      -1.9416    -2.8622     -3.413
  %  1% C.V       -2.5672    -3.4329    -3.9623
  %  P-value      0.033353      0.001      0.001


%KPSS test..Similar to the pptest
[h4,pValue4,stat4,cValue4,reg4]=kpsstest(data, 'alpha', alpha,'lags', ...
    Bandwidth,  'trend', false);
[h5,pValue5,stat5,cValue5,reg5]=kpsstest(data, 'alpha', alpha,'lags', ...
    Bandwidth, 'trend', true);
%Prepare a table to present the test outputs
kpsstable=[stat4(1) cValue4 pValue4(1); stat5(1) cValue5 pValue5(1)];
kpsstable
%  0.9964    0.3470    0.4630    0.7390    0.0100
%  1.0041    0.1190    0.1460    0.2160    0.0100
kpssout=array2table(kpsstable','RowNames', {'Test Stat' ' 10% C.V.' '5% C.V.' '1% C.V.' 'P-value' },...
'VariableNames', {'Drift' 'Trend'}); 
%
%               Drift     Trend 
%                 _______    ______

%    Test Stat    0.99641    1.0041
%    10 % C.V.      0.347     0.119
%    5% C.V.        0.463     0.146
%    1% C.V.        0.739     0.216
%    P-value         0.01      0.01

%Print the results in the command window
disp(strcat('ADF test with optimal lag selection using', {' '}, sel_method));
% ADF test with optimal lag selection using BIC'
disp(dfout);
  %                 None       Drift      Trend 
  %               ________    _______    _______
 %
 %   Test Stat     -2.0757    -5.3246    -5.3251
 %   10% C.V.      -1.6173    -2.5672    -3.1279
 %   5% C.V.       -1.9416    -2.8622     -3.413
 %   1% C.V.       -2.5672    -3.4329    -3.9623
 %   P-value      0.036608      0.001      0.001
%    Lag used           10         10         10
disp('PP test with Newey-West fixed Bandwidth');
% PP test with Newey-West fixed Bandwidth
disp(ppout);
%                   None       Drift      Trend 
%                 ________    _______    _______

%    Test Stat     -2.1144    -6.6476    -6.6493
%    10% C.V.      -1.6173    -2.5672    -3.1279
%    5% C.V.       -1.9416    -2.8622     -3.413
%    1% C.V.       -2.5672    -3.4329    -3.9623
%    P-value      0.033353      0.001      0.001
disp('KPSS test with Newey-West fixed Bandwidth');
% KPSS test with Newey-West fixed Bandwidth
disp(kpssout);
%                  Drift     Trend 
%                 _______    ______
%
%    Test Stat    0.99641    1.0041
%    10% C.V.       0.347     0.119
%    5% C.V.        0.463     0.146
%    1% C.V.        0.739     0.216
%    P-value         0.01      0.01

