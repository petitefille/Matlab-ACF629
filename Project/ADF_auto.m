function out=ADF_auto(data,dif,max_lags,sel_method)
% dfout=ADF_auto(retOtC,0,30,sel_method);
%This function computes the augmented Dickey Fuller (ADF test with automatic lag
%selection based on an information criteria from AIC, BIC and HQIC. Three
%sets of tests are performed: the ADF test with no drift or trend, with a
%drift and with both drift and trend. 

%Inputs: 
%data: T-by-1 matrix of data => must be a vector


%dif: non-negative integer, the ADF test is performed after differencing
%(0)
%the original series 'dif' times. (30)
%max_lags: non-negative integer, the maximum lags included in the ADF test
%sel_method = selection_method: choose from 'AIC', 'BIC' and 'HQIC', the information
%criteria minimized to select an optimal lag length

%output: A 6-by-3 table documenting 
% the test statistics
% the critical values (corresponding to respective significance level)
% p-value and
%the optimal lags chosen


%Checking selection methods
sel_method_choice={'AIC' 'BIC' 'HQIC'};

% if sel_method not i sel_method_choice (ismember() == 0)
if ismember(sel_method,sel_method_choice)==0 
    % ismember(x,y)
    % returnerer 1 hvis el x er i liste y
    % 0 ellers
    error('Selection method not supported.');
end
%Checking input type of lags (must be integers)
% rem(a,b) returns remainder after division of a by b (integer division)
% 0.1/1 = 0*1 + 0.1 so rem(0.1,1) = 0.1
% ~= is not equal to
if rem(dif,1)~=0 || rem(max_lags,1)~=0
      error('Wrong input type.');
end
%Differencing the data
sdata=data; % sdata er nå en kopi av data
while dif>0 % så lenge diff er større enn 0
sdata=diff(data,dif); % sdata is always equal to 1st difference ?
dif=dif-1;
end

warning('off');
%Do an ADF test for 0:max_lags with no drift or trend 
[h,pValue,cValue,stat, regar]=adftest(sdata,'Model','AR', 'lags', 0:max_lags);
% Y = sdata = univariate time series
% 'alpha' significance level 0.05
% 'lags': nr of lagged difference terms to include in the model

%Do an ADF test for 0:max_lags with drift but no trend 
% here I will obtained values for each lag
[hd,pValued,cValued,statd,regard] = adftest(sdata,'Model','ARD','lags',0:max_lags);
% hd: 1 x 31
% pValued: 1 x 31
% cValued: 1 x 31 
% statd: 1 x 31
% regard: see matlab (multi dimensional)
%Do an ADF test for 0:max_lags with drift and trend
[hts,pValuedts,cValuedts,statts,regts] = adftest(sdata,'Model','TS','lags',0:max_lags);
%Store the information criteria in three different matrices
ICAR=[(regar.AIC); (regar.BIC); (regar.HQC)]; % 3 x 31, row = selection method
ICARD=[(regard.AIC); (regard.BIC); (regard.HQC)];
ICTS=[(regts.AIC); (regts.BIC); (regts.HQC)];
% Choosing the minimum information criteria from the corresponding matrix
[~, minICARi]=min(ICAR,[],2);
              % is a 3 x 1 vector containing the index of the
              % min value of each row (3
              % rows) with respect to which method is used
              % AIC, BIC or HQC
[~, minICARDi]=min(ICARD,[],2);
[~,minICTSi]=min(ICTS,[],2);
%Initialize the matrix Qtable to store results
Qtable=zeros(3,6);
%Choose the corresponding best lag according to the selected information
%criterion

% tf = strcmp(s1,s2) compares s1 and s2 and returns 1 (true) if the two are identical and 0 (false) otherwise. 
% Text is considered identical if the size and content of each are the same. The return result tf is of data type logical.
if strcmp(sel_method,'AIC')==1 % sel_method = 'AIC'
    lagar=minICARi(1);
    lagard=minICARDi(1);
    lagts=minICTSi(1);
elseif  strcmp(sel_method,'BIC')==1
    lagar=minICARi(2); % 11
    lagard=minICARDi(2); % 11
    lagts=minICTSi(2); % 11
else
    lagar=minICARi(3);
    lagard=minICARDi(3);
    lagts=minICTSi(3);
end
% Use the chosen lag to produce ADF test results
[h1,p1,df1,cv1,reg1] = adftest(sdata,'model','AR','alpha',[0.1 0.05 0.01],'lags',lagar);
[h2,p2,df2,cv2,reg2] = adftest(sdata,'model','ARD','alpha',[0.1 0.05 0.01],'lags',lagard);
[h3,p3,df3,cv3,reg3] = adftest(sdata,'model','TS','alpha',[0.1 0.05 0.01],'lags',lagts);
%Store the relevant statistics in a table
% 1st row, 1st col
Qtable(1,1)=df1(1); % test statistic as standard test statistic
% 1st row and 2 col
Qtable(2,1)=df2(1);
Qtable(3,1)=df3(1);
% 1st column, row 2-4 
Qtable(1,2:4)=cv1; % critical values 
% 2nd col, row 2-4
Qtable(2,2:4)=cv2; % critical values 
% 3rd col, row 2-4 
Qtable(3,2:4)=cv3; % critical values
% 1st col, 5th row 
Qtable(1,5)=p1(1); % p value with respect to standard test statistic
% 2nd col, 5th row
Qtable(2,5)=p2(1);
% 3rd col, 5th row
Qtable(3,5)=p3(1);
% 1st col, 6th row
Qtable(1,6)=lagar-1; % appropriate nr of lags 
% 2nd col, 6th row
Qtable(2,6)=lagard-1;
% 3rd col, 6th row
Qtable(3,6)=lagts-1;
%Prepare the output.
% Here: Qtable is transposed!!!
out=array2table(Qtable','RowNames', {'Test Stat' ' 10% C.V.' '5% C.V.' '1% C.V.' 'P-value' 'Lag used'},...
'VariableNames', {'None' 'Drift' 'Trend'}); 
end
