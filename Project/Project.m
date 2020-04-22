% Coursework

% Part 1: ARMA-GARCH Model
% a) For your chosen stock, compute the tick-by-by tick log-return and the
% open-to-close log return series for each trading day in the dataset.

clearvars; clc;

% set working directory 
% cd('C:\Users\Eier\Documents\Lancaster\ACF609\Project')

%load csv file into MATLAB as table
dmat=readtable('IBM_trade_quote.csv');

% extract data
data=dmat{:,1:4};

size(dmat); % 10535140           6

%extract all dates 
dateTbT = data(:,2); % 10535140 x 1 vector

% convert dates to str
dateTbT_str = num2str(dateTbT);

% convert str to datetime
dateTbT_dt = datetime(dateTbT_str,'InputFormat','yyyyMMdd','Format','dd.MM.yyyy');

%Extract unique day counters
dc=unique(data(:,1)); % 2517 x 1

% Extract unique dates
ddate=unique(data(:,2));
% ddate: 2517 x 1 matrix 

% convert unique dates to string
dateOtC_str = num2str(ddate);

% convert unique dates from string to datetime format
dateOtC_dt = datetime(dateOtC_str,'InputFormat','yyyyMMdd','Format','dd.MM.yyyy');


% initializing the open-to-close log return vector
retOtC=zeros(length(dc),1);

% initializing the tick-by-tick log return vector
retTbT=diff(log(data(:,4))); 

%For loop for all trading days
for i=1:length(dc) 
   pricet=data(data(:,1)==dc(i),3:4);
     
   % collecting open-to-close log returns
   retOtC(i)=log(pricet(end,2))-log(pricet(1,2));
   
end

% length(retTbT) % 10535139 x 1
% length(retOtC) % 2517 x 1 

% b) Provide a descriptive analysis for both return series computed above.
% Discuss your findings and relate them to the stylized facts of financial
% data.

%get the size of the data
[rowsTbT, colsTbT] = size(retTbT); % [ 10535139, 1 ]
[rowsOtC, colsOtC]=size(retOtC); % [2517,1] 

names= {'TickByTick','OpenToClose'};

%Computing the statistics

%mean
meandata= [mean(retTbT), mean(retOtC)];

%variance
vardata= [var(retTbT), var(retOtC)];

% minimum 
mindata= [min(retTbT), min(retOtC)];

% maximum
maxdata = [max(retTbT),max(retOtC) ];

%standard deviation.
stddata= [std(retTbT), std(retOtC) ];
   
% skewness = the sample skewness of X
skedata= [ skewness(retTbT), skewness(retOtC)];

%kurtosis
% sample kurtosis 
kurdata=[kurtosis(retTbT), kurtosis(retOtC) ];

%quantiles
qvector=[.01 .05 .1 .5 .9 .95 .99]; % vektor 
qdata=[quantile(retTbT,qvector); quantile(retOtC,qvector)];

%Prepare a table for the results
%combine the statistics first:
stat=[meandata; vardata; maxdata;mindata; stddata; skedata; kurdata;qdata'];
    
%convert it to table with proper titles
stattable=array2table(stat,'VariableNames',names,'RowNames', {'Mean' ...
    'Variance' 'Maximum' 'Minimum' 'Std.Dev.' 'Skewness' 'Kurtosis' ...
    'Q(0.01)' 'Q(0.05)' 'Q(0.1)' 'Q(0.5)' 'Q(0.9)' 'Q(0.95)' 'Q(0.99)'});

%display the table to the output window
disp(stattable);
%                TickByTick     OpenToClose
%                ___________    ___________
%
%    Mean         4.5862e-08    0.00097579 
%    Variance     6.0638e-08    0.00012392 
%    Maximum        0.074876      0.066563 
%    Minimum       -0.079237     -0.066163 
%    Std.Dev.     0.00024625      0.011132 
%    Skewness        -7.4979      -0.24509 
%    Kurtosis         6533.7        7.7745 
%    Q(0.01)     -0.00060646     -0.031155 
%    Q(0.05)     -0.00029193     -0.016138 
%    Q(0.1)      -0.00019445     -0.010567 
%    Q(0.5)                0    0.00086823 
%    Q(0.9)       0.00019389      0.012998 
%    Q(0.95)      0.00029468      0.017482 
%    Q(0.99)      0.00061121      0.029684 


% line graphs of returns 
plotsize=[100,100,1200,550]; 
fig1=figure;
subplot(1,2,1);
plot(dateTbT_dt(2:end),diff(log(data(:,4)))); 
title('Tick-by-tick return');   
subplot(1,2,2);
plot(dateOtC_dt,retOtC);
title('Open-to-close return');
suptitle('Line graphs of IBM returns');

% histograms
fig1.Position=plotsize; 
fig2=figure;
subplot(1,2,1);
histogram(retTbT);
title('Tick-by-tick return');
subplot(1,2,2);
histfit(retOtC);
title('Open-to-close return');
suptitle('Histograms of IBM returns');
fig2.Position=plotsize;

%Correlogram and partial correlogram up to lag L 
L=36;
fig3=figure;
Q=zeros(L,2);
[acfTbT,lag1TbT,bound1TbT]=autocorr(retTbT,L);   
[pacfTbT,lag2TbT,bound2TbT]=parcorr(retTbT,L); 
%compute Ljung-Box test statistics 
Q(:,1)=length(retTbT)*(length(retTbT)+2)*cumsum((acfTbT(2:L+1).^2)./(length(retTbT)-lag1TbT(2:L+1)));
subplot(1,2,1);    
hold on;   
p1=stem(1:1:L,acfTbT(2:L+1),'filled'); %ploting acf
p2=stem(1:1:L,pacfTbT(2:L+1)); %ploting pacf
p=plot(1:1:L,[ones(L,1).*bound1TbT(1) ones(L,1).*bound1TbT(2)]); %ploting confidence bounds
hold off; 
%Graph settings
p1.Color='r';
p1.MarkerSize=4;
p2.Color='b';
p2.MarkerSize=4;
p(1).LineWidth=1;
p(1).LineStyle='--';
p(1).Color='k';
p(2).LineWidth=1;
p(2).LineStyle='--';
p(2).Color='k';
title('Tick-by-tick return'); 

[acfOtC,lag1OtC,bound1OtC]=autocorr(retOtC,L);   
[pacfOtC,lag2OtC,bound2OtC]=parcorr(retOtC,L);
%compute Ljung-Box test statistics 
Q(:,2)=length(retOtC)*(length(retOtC)+2)*cumsum((acfOtC(2:L+1).^2)./(length(retOtC)-lag1OtC(2:L+1)));   
subplot(1,2,2);    
hold on;  
p1=stem(1:1:L,acfOtC(2:L+1),'filled'); %ploting acf
p2=stem(1:1:L,pacfOtC(2:L+1)); %ploting pacf
p=plot(1:1:L,[ones(L,1).*bound1OtC(1) ones(L,1).*bound1OtC(2)]); % ploting confidence bounds
hold off; 
%Graph settings
p1.Color='r';
p1.MarkerSize=4;
p2.Color='b';
p2.MarkerSize=4;
p(1).LineWidth=1;
p(1).LineStyle='--';
p(1).Color='k';
p(2).LineWidth=1;
p(2).LineStyle='--';
p(2).Color='k';
title('Open-to-close return');     
%legend settings
legend('ACF', 'PACF','95% Confidence Bounds', 'Location','southoutside','Orientation','horizontal'); 
suptitle('Correlogram and partial correlogram of IBM returns');
fig3.Position=plotsize;

% ACF for IBM squared returns 
L=36;
fig4=figure;
[acfTbT,lag1TbT,bound1TbT]=autocorr(retTbT.^2,L);
subplot(1,2,1);    
hold on; 
p1=stem(1:1:L,acfTbT(2:L+1),'filled'); %ploting acf
p=plot(1:1:L,ones(L,1).*bound1TbT(1)); %ploting confidence bound
hold off; 
%Graph settings
p1.Color='r';
p1.MarkerSize=4;
p.LineWidth=1;
p.LineStyle='--';
p.Color='k';
title('tick-by-tick squared returns'); 

[acfOtC,lag1OtC,bound1OtC]=autocorr(retOtC.^2,L);   
subplot(1,2,2);    
hold on;  
p1=stem(1:1:L,acfOtC(2:L+1),'filled'); 
p=plot(1:1:L,ones(L,1).*bound1OtC(1)); 
hold off; 
%Graph settings
p1.Color='r';
p1.MarkerSize=4;
p.LineWidth=1;
p.LineStyle='--';
p.Color='k';
title('Open-to-close squared returns');     
%legend settings
legend('ACF','95% Upper Confidence Bound', 'Location','southoutside','Orientation','horizontal');  
suptitle('ACF of IBM squared returns');
fig4.Position=plotsize;

% ----------------------------------------------------------
%Jacque-Bera test for normality
% H0: Normal
JB = [(length(retTbT))/6*(skedata(1).^2+0.25*(kurdata(1)-3).^2), (length(retOtC))/6*(skedata(2).^2+0.25*(kurdata(2)-3).^2)]; % data: log returns
JBpvalue=chi2cdf(JB,2,'upper'); % p-verdi

%Ljung-Box test for correlation (pvalue computed using the Q obtained from the
%correlogram)
LBpvalue=chi2cdf(Q,repmat((1:1:L)',1,2),'upper'); 
teststat=[ JB;Q]; 
pstat=[ JBpvalue; LBpvalue]; 
testtable=cell(L+1,3); 
%Use * to indicate significance: ***: significant at 1%, **: significant at
%5%, *: significant at 10%


for i=1:L+1 
    if i==1 
        testtable(i,1)={'JB stat'}; 
    else
        testtable(i,1)={strcat('LB(', num2str(i-1,'%5.0f'),')')};
            
    end
    for j=2:3   
        %  p value <= alpha: reject H0
        if pstat(i,j-1)<=0.01  % H0 blir forkastet hvis alpha = 0.01 
            testtable(i,j)={strcat(num2str(teststat(i,j-1),'%5.2f'),'***')};% test verdi og ikke p value
                                                                            % settes inn 
        elseif pstat(i,j-1)<=0.05
              testtable(i,j)={strcat(num2str(teststat(i,j-1),'%5.2f'),'**')};
        elseif pstat(i,j-1)<=0.1
              testtable(i,j)={strcat(num2str(teststat(i,j-1),'%5.2f'),'*')};
        else
              testtable(i,j)={num2str(teststat(i,j-1),'%5.2f')}; % H0 blir ikke forkastet 
                                                                 % ingen
                                                                 % stjerner
                                                                 % 
        end
    end
end
%create table of results
testtable=cell2table(testtable(:,2:end), 'VariableNames',  names, 'RowNames',testtable(:,1));

disp(testtable); 

%                     TickByTick          OpenToClose 
%               ______________________    ____________
%
%    JB stat    '18722109883869.32***'    '2415.96***'
%    LB(1)      '25761.20***'             '2.65'      
%    LB(2)      '25778.43***'             '3.11'      
%    LB(3)      '25779.20***'             '16.09***'  
%    LB(4)      '25816.37***'             '18.43***'  
%    LB(5)      '25818.34***'             '18.43***'  
%    LB(6)      '25833.15***'             '23.36***'  
%    LB(7)      '25838.26***'             '23.62***'  
%    LB(8)      '25850.35***'             '25.50***'  
%    LB(9)      '25869.65***'             '30.06***'  
%    LB(10)     '25886.20***'             '33.37***'  
%    LB(11)     '25888.52***'             '34.39***'  
%    LB(12)     '25897.17***'             '34.93***'  
%    LB(13)     '25897.37***'             '35.01***'  
%    LB(14)     '25900.00***'             '51.46***'  
%    LB(15)     '25902.92***'             '52.22***'  
%    LB(16)     '25916.81***'             '53.76***'  
%    LB(17)     '25967.60***'             '54.41***'  
%    LB(18)     '26009.23***'             '56.17***'  
%    LB(19)     '26009.44***'             '56.23***'  
%    LB(20)     '26012.37***'             '58.16***'  
%    LB(21)     '26013.39***'             '58.64***'  
%    LB(22)     '26022.32***'             '61.00***'  
%    LB(23)     '26022.41***'             '63.76***'  
%    LB(24)     '26028.16***'             '63.99***'  
%    LB(25)     '26031.30***'             '64.85***'  
%    LB(26)     '26050.02***'             '70.09***'  
%    LB(27)     '26091.34***'             '71.05***'  
%    LB(28)     '26092.77***'             '71.17***'  
%    LB(29)     '26093.54***'             '72.00***'  
%    LB(30)     '26118.14***'             '72.36***'  
%    LB(31)     '26119.09***'             '74.24***'  
%    LB(32)     '26119.44***'             '75.36***'  
%    LB(33)     '26140.17***'             '75.74***'  
%    LB(34)     '26203.04***'             '76.08***'  
%    LB(35)     '26211.03***'             '88.51***'  
%    LB(36)     '26232.15***'             '92.00***'


% c) For the daily open-to-close return time series computed in a), find
% the best model in the general ARMA-GARCH modelling framework. You should
% provide  detailed model selection procedures and justify your choice 
% Use a sensible set of competing ARMA-GARCH specifications.

%Set the significance level of the test
alpha=[0.1,0.05,0.01];

%Augmented Dickey-Fuller test with automatic lag selection based on
%information criterion.

sel_method='BIC';

dfout=ADF_auto(retOtC,0,30,sel_method);% dif = 0
     
%For PP and KPSS test, the selection of lags is more involved. Here I used
%a fixed bandwith for the Newey-West type variance-covariance estimates,
%which may not be the best selection method for the two type of tests. 
Bandwidth=floor(4* (length(retOtC)^(2/9)) ); % 22 

%Philip Perron test
%Perform three sets of tests for three different model specifications
[h1,pValue1,stat1,cValue1,reg1]=pptest(retOtC,'Model','AR', 'alpha', alpha,'lags', Bandwidth);
[h2,pValue2,stat2,cValue2,reg2]=pptest(retOtC,'Model','ARD', 'alpha',alpha, 'lags', Bandwidth);
[h3,pValue3,stat3,cValue3,reg3]=pptest(retOtC,'Model','TS',  'alpha', alpha,'lags', Bandwidth);

%Prepare a table to present the test outputs
pptable=[stat1(1) cValue1 pValue1(1); stat2(1) cValue2 pValue2(1) ;stat3(1) cValue3 pValue3(1)];
ppout=array2table(pptable','RowNames', {'Test Stat' ' 10% C.V.' '5% C.V.' '1% C.V.' 'P-value' },...
'VariableNames', {'None' 'Drift' 'Trend'}); 

%Print the results in the command window
disp(strcat('ADF test with optimal lag selection using', {' '}, sel_method));
% ADF test with optimal lag selection using BIC'
disp(dfout);

 %                 None       Drift      Trend 
 %                _______    _______    _______
%
%    Test Stat     -34.09    -34.463    -34.467
%    10% C.V.     -1.6168    -2.5683    -3.1283
%    5% C.V.      -1.9416    -2.8642    -3.4141
%    1% C.V.      -2.5689    -3.4366    -3.9675
%    P-value        0.001      0.001      0.001
%    Lag used           0          0          0
    
disp('PP test with Newey-West fixed Bandwidth');
% PP test with Newey-West fixed Bandwidth
disp(ppout);
%                   None       Drift      Trend 
%                _______    _______    _______
%
%    Test Stat    -49.839    -49.608    -49.593
%    10% C.V.     -1.6168    -2.5683    -3.1283
%    5% C.V.      -1.9416    -2.8642    -3.4141
%    1% C.V.      -2.5689    -3.4366    -3.9675
%    P-value        0.001      0.001      0.001

T=length(retOtC); % 2517

%Picking the best model from ARMA(0,0) to ARMA(pmax,qmax) based on information
%criteria

%initiallizing information criteria mat
pmax=1;
qmax=1;
AICmat=zeros(pmax+1,qmax+1); 
BICmat=zeros(pmax+1,qmax+1);
HQICmat=zeros(pmax+1,qmax+1);

%Compute the log-likelihood and information criteria for the models considered
%Note that I started with ARMA(0,0) and loop over all the models to
%ARMA(pmax,qmax).
for i=1:pmax+1 % pmax = 0,1
    for j=1:qmax+1 % qmax = 0,1
        Mdltemp=arima(i-1,0,j-1);
        [~,~,logL,~] = estimate(Mdltemp,retOtC,'Display','off');
             
        %You can also use the aicbic function to compute AIC or BIC.
        AICmat(i,j)=-2*logL/T+2*(-1+i+j)/T; 
        % i+j-1 
        BICmat(i,j)=-2*logL/T+(-1+i+j)*log(T)/T;
        HQICmat(i,j)=-2*logL/T+2*(-1+i+j)*log(log(T))/T;
    end
end

%Find the location of the minimum information critieria. 
[paic,qaic]=find(AICmat==min(min(AICmat)));  
                                            
[pbic,qbic]=find(BICmat==min(min(BICmat)));
[phqic,qhqic]=find(HQICmat==min(min(HQICmat)));

%Choose best model using BIC for example
pbest=pbic-1; % 1
qbest=qbic-1; % 1

% [P,Q] = ARMA_optimal(retOtC,pmax,qmax,'BIC');% P = 1, Q = 1

Mdlbest=arima(pbest,0,qbest);
[EstMdlbest,EstParamCov,logL1,info] = estimate(Mdlbest,retOtC);
%EstMdlbest: model containing parameter estimates
 
% ARIMA(1,0,1) Model (Gaussian Distribution):
 
%              Value       StandardError    TStatistic      PValue   
%            __________    _____________    __________    ___________
%  Constant    8.0012e-05      3.267e-05        2.4491         0.014322
%  AR{1}          0.92247       0.024813        37.177               0
%  MA{1}          -0.8914       0.030111       -29.603               0
% Variance    0.00012286     1.9245e-06        63.841                0

%The parameter estimates are stored in the structure EstMdlbest. To retrieve
%them, use the following command:
bhat=[EstMdlbest.Constant; EstMdlbest.AR{1};  EstMdlbest.MA{1}; EstMdlbest.Variance   ];
  
% The standard errors can be retrieved from the EstParamCov matrix
se=sqrt(diag(EstParamCov)); 
 
%construct t-statistics and p-values ourselves:
tstat=bhat./se; % tstat = estimated value of parameters/ SE(parameters)

pvalue=(1-normcdf(abs(tstat)))*2;  
%Prepare a table for the estimation results ourselves. 
armadata=[bhat se tstat pvalue ];
armaout=array2table(armadata, 'RowNames', {'c' 'phi' 'theta' 'sigmasq'}, 'VariableNames', {'bhat' 'SE' 'tstat' 'pvalue'});
disp(armaout);
%                  bhat           SE         tstat      pvalue 
%               __________    __________    _______    ________
%
%    c          8.0012e-05     3.267e-05     2.4491    0.014322
%    phi           0.92247      0.024813     37.177           0
%    theta         -0.8914      0.030111    -29.603           0
%    sigmasq    0.00012286    1.9245e-06     63.841           0

% Residual analysis 
[E,V,logL2]=infer(EstMdlbest,retOtC);

%ARCH LM test to check for potential heteroscedasticity in the  error
%term. We only check for the first lag.
 [h,pValue,stat,cValue] = archtest(E,'Lags',1);
  % h = 1
  % pValue = 0


%Compute optimal GARCH lags using a function GARCH_optimal.m. 
[gbest,abest,ICmat] = GARCH_optimal(retOtC, 2, 2,1, 1, 'BIC');
% gbest = 1
% abest = 1

%Using the optimal GARCH lags 
Mdl=arima('ARLags',1','MALags',1,'Variance',garch(gbest,abest));

%Estimate the model
[estmdl, estParamCov,logL, info] =estimate(Mdl,retOtC);
 
% ARIMA(1,0,1) Model (Gaussian Distribution):
 
%                  Value       StandardError    TStatistic      PValue   
%                __________    _____________    __________    ___________
%
%    Constant    0.00014307     5.6657e-05        2.5252         0.011563
%    AR{1}          0.88452        0.04154        21.293                0
%    MA{1}         -0.84785        0.04756       -17.827                0

 
 
%    GARCH(1,1) Conditional Variance Model (Gaussian Distribution):
 
%                  Value       StandardError    TStatistic      PValue  
%                __________    _____________    __________    __________

%    Constant    3.0732e-06     7.9474e-07         3.867      0.00011019
%    GARCH{1}       0.87601       0.012328        71.058               0
%    ARCH{1}       0.095708       0.009778        9.7881               0


% d) Perform diagnostic checks on your best ARMA-GARCH model. Discuss the
% goodness-of-fit of your model.


% residuals for the ARMA(1,1)-GARCH(1,1) fitted model
[E,V,logl]=infer(estmdl,retOtC);

%Compute fitted standardized error term
z=E./sqrt(V);

%Diagnostic testing
%testing against a standard normal distribution
qqplot(z);
% Jarque Bera test
[h,p,stat,cv]=jbtest(z); % testing standardized residuals
% h: 1 => not normally distributed
% p: 1.0000e-03

%Detecting remaining ARCH effects
[h2,p2,stat2,cv2] = archtest(z,'Lags',1:10);
% h2: 0   0   0   0   0   0   0   0   0   0 0 
% p:   0     0     0     0     0     0     0     0     0     0
% no conditional heteroscedasticity
% our model captures the heteroscedasticity present in data


%t-distributed GARCH models
%Changing our definition of Mdl slightly
Mdl2=Mdl;% Mdl.Distribution = "Gaussian"
Mdl2.Distribution='t';
%Estimate the model
[estmdl2, EstParamCov2, logL2, info2]=estimate(Mdl2,retOtC);

% residuals for ARMA(1,1)-GARCH(1,1) model (t-distributed)
[E2,V2,logl2]=infer(estmdl2,retOtC);

%Compute fitted standardized error term
z2=E2./sqrt(V2);
%We need to make sure that our specification using the t-distribution is 
%correct. We can use a Q-Q plot to check it.
%Firstly, we get the estimated DoF and define a t-distribution with the
%estimated DoF
dof=estmdl2.Distribution.DoF; % 8.3607
dist=makedist('tLocationScale',0,1,dof);

%We can then plot the standardized residual accordingly
% to compare the distribution of our standardized residuals and
% the t distribution 
qqplot(z2,dist);
%It is clear that the goodness-of-fit is much better compared to a normal
%assumption.

%Asymmetric GARCH. We consider only EGARCH
mdl3=arima('ARLags',1,'MALags',1,'Variance',egarch(1,1));
[estmdl3, EstParamCov3,logL3, info3]=estimate(mdl3,retOtC); 

[E3,V3,logL3]=infer(estmdl3,ret); % residuals EGARHC(1,1)

%Compute fitted standardized error term
z3=E3./sqrt(V3);
% assessing EGARCH model 
qqplot(z3);

%GJR-GARCH:
mdl4=arima('ARLags',1,'MALags',1,'Variance',gjr(1,1)); 
[estmdl4,EstParamCov4,logL4,info4] =estimate(mdl4,retOtC); 

[E4,V4,logl4]=infer(estmdl4,retOtC);% residuals 
z4=E4./sqrt(V4); % standardized residuals
qqplot(z4);

% GJR-GARCH specified with a t-distirbuted error term.
mdl6 = mdl4; 
mdl6.Distribution='t';
[estmdl6, EstParamCov6, logL6, info6]=estimate(mdl6,retOtC); 

 
 %   ARIMA(1,0,1) Model (t Distribution):
 
 %                 Value       StandardError    TStatistic      PValue  
 %               __________    _____________    __________    __________
%
%    Constant    0.00013195      5.956e-05        2.2155        0.026729
%    AR{1}          0.87059       0.051299        16.971               0
%    MA{1}         -0.84175       0.057336       -14.681               0
 %   DoF             9.1731         1.6395         5.595      2.2063e-08

 
 
%    GJR(1,1) Conditional Variance Model (t Distribution):
 
 %                    Value       StandardError    TStatistic      PValue  
 %                  __________    _____________    __________    __________

 %   Constant       2.5798e-06     8.5835e-07        3.0055       0.0026511
 %   GARCH{1}          0.89248       0.014831        60.178               0
 %   ARCH{1}          0.032586       0.013487         2.416        0.015691
  %  Leverage{1}      0.099503       0.020197        4.9267      8.3639e-07
  %  DoF                9.1731         1.6395         5.595      2.2063e-08

[E6,V6,logl6]=infer(estmdl6,retOtC);
z6=E6./sqrt(V6); % standardized residuals

dof=estmdl6.Distribution.DoF;
dist=makedist('tLocationScale',0,1,dof);
qqplot(z6,dist);

%Since the GJR-GARCH model nests the GARCH model, we can conduct a
%likelihood-ratio test to see if the extra parameter in the GJR-GARCH model
%improves the likelihood significantly.
% first: compare ARMA(1,1)-GARCH(1,1) and ARMA(1,1)-GJR(1,1) 
 % [h3,p3,stat3,cv3] = lratiotest(logL4,logl,1);
 % logl4:  unrestricted model
 % logl: restricted model
 % h3 = 1; p3 = 1.8839e-09 
 % h3 = 1 means unrestricted model (GJR-GARCH) is better
%The test is significant at any conventional significance level, suggesting
%that we should use an asymmetric structure when modelling financial
%returns with GARCH-type structure.

% 2nd: likelihood ratio test comparing ARMA(1,1)-GARCH(1,1) and ARMA(1,1)_GJR(1,1) (t-distributed)
[h4,p4,stat4,c4] = lratiotest(logl6,logl,1);
% h4 : 1
% p4:  0

% Jarque Bera test for the ARMA(1,1)-GJR(1,1) (t-distributed)
[h7,p7,stat7,cv7]=jbtest(z6); % testing standardized residuals
% h7: 1 => not normally distributed
% p7: 1.0000e-03

%Detecting remaining ARCH effects for the ARMA(1,1)-GJR(1,1) (t-distributed)
[h8,p8,stat8,cv8] = archtest(z6,'Lags',1:10);
% h8: 0   0   0   0   0   0   0   0   0   0 0 => no conditional
% heteroscedasticity
% our model captures the heteroscedasticity present in data

% Residual analysis for for the ARMA(1,1)-GJR(1,1) (t-distributed)
plotsize=[100,100,1200,550];
fig1 = figure;
hold on;
lag=100;
[acf1,~,bound]=autocorr(E6,lag); % 101 x 1 (from la 0 to lag 100)
[acf2,~,bound1]=autocorr(E6.^2,lag); % bound and bound1 are the same
p1=stem(1:1:lag,acf1(2:lag+1),'filled'); % ploting acf in blue and filled
p2=stem(1:1:lag,acf2(2:lag+1)); %ploting acf^2 not filled in red
p=plot(1:1:lag,[ones(lag,1).*bound1(1) ones(lag,1).*bound1(2)]); %ploting confidence bounds
hold off;
%Graph settings
p1.Color='r';
p1.MarkerSize=4;
p2.Color='b';
p2.MarkerSize=4; 
p(1).LineWidth=2;
p(1).LineStyle=':'; 
p(1).Color='k'; 
p(2).LineWidth=2;
p(2).LineStyle=':';
p(2).Color='k'; 
title('Correlogram of the residual and squared residual');
legend('Residual', 'Squared residual','95% Confidence Bounds', 'Location','southoutside','Orientation','horizontal');


% Part 2: Realized Volatility signature plot

% e) Construct a volatility signature plot with daily annualized realized
% volatility estimators based on the high-frequency data for the chosen
% stock.
% cd('C:\Users\Eier\Documents\Lancaster\ACF609\Lab6')
% right click on mat folder in Lab6 (ACF609) add to path (folders + subfolders)

%Set the time grid to compute the volatility signature plot
 tgrid=1:5:1800; % plot volatility signature plot = sampling intervals in sec

 %Initialize a matrix for  volalatility signature plot
 VSP=zeros(length(dc),length(tgrid)); % store volatility signature plo
 
 % 252 is the nr of trading days in a year
 
 %initialize various measures of RV
 rv=zeros(length(dc),1); % 252 x 1 matrise med elementer 0 
 

%For loop for all trading days
for i=1:length(dc) 
   pricet=data(data(:,1)==dc(i),3:4);
   
%Using the MFEtoolbox to construct volatility signature plot based on the
%sampling frequencies in the tgrid
for j=1:length(tgrid) 
   RV=realized_variance(pricet(:,2),pricet(:,1),...
       'seconds','CalendarTime',tgrid(j)); 
   VSP(i,j)=RV*252;
end
end

% f) Select an optimal sampling frequency and construct the daily
% annualized realized volatility estimator based on the optimal sampling
% frequency. 
 
%For loop for all trading days
for i=1:length(dc)  
   pricet=data(data(:,1)==dc(i),3:4);
   
%Compute the 700-sec sample RV and
[RV,~]=realized_variance(pricet(:,2),pricet(:,1),...
       'seconds','CalendarTime',700,5); % split into 5 groups (k: 60)
                                        % 5 is nr of subsamples
                                        % 300 / 60 = 5 min !
                                        % every 300 min = 5 min
rv(i)=RV*252;%700 sec simple RV for day i (annualized)

end

figure('units','normalized','outerposition',[0 0 1 1]);
%Plot the volatility signature plot
plot(tgrid,mean(VSP)); % VSP
title('Volatility Signature Plot');
xlabel('Sampling Interval In Seconds');
ylabel('Average RV');

figure('units','normalized','outerposition',[0 0 1 1]);
%Plot RV measure
plot(rv);
legend({'700 sec RV'});

% g) Split the dataset into an in-sample period(03/01/2005 to 31/12/2012)
% and an out-of-sample period (02/01/2013 to 31/12/2014). Compute the
% rolling-window , one-step-ahead (one-day-ahead) volatility forecasts with
% the ARMA-GARCH model and a HAR-RV specification. To reduce computational
% burden, estimate the ARMA-GARCH and HAR-RV model parameters only once on
% the basis of the initial in-sample period. 

T=length(retOtC);
%Set an insample period 
insmpl=data(year(dateTbT_dt)<=2012,1:4); 
%Similarly, set an out-of-sample period
outofsmpl=data(year(dateTbT_dt)>2012,1:4); 

% compute open to close log returns for insample and out of sample
dcIn=unique(insmpl(:,1));
dcOut = unique(outofsmpl(:,1));

retOtCIn=zeros(length(dcIn),1);
retOtCOut=zeros(length(dcOut),1);

for i=1:length(dcIn) 
   pricet=insmpl(insmpl(:,1)==dcIn(i),3:4);
   retOtCIn(i)=log(pricet(end,2))-log(pricet(1,2));
   
end

for i=1:length(dcOut) 
   pricet=outofsmpl(outofsmpl(:,1)==dcOut(i),3:4);
   retOtCOut(i)=log(pricet(end,2))-log(pricet(1,2)); 
end

% cd('C:\Users\Eier\Documents\Lancaster\ACF609\Lab4');

% optimal model: ARMA(1,1)-GJR(1,1) (t-distributed)
mdl=arima('ARLags',1,'MALags',1,'Variance',gjr(1,1));
mdl2 = mdl; 
mdl2.Distribution='t';
[estmdl, EstParamCov, logL, info]=estimate(mdl2,retOtC);
%Get the names of the series
RVnames = {'RV','GJR-GARCH-t'};

% [y, ymse, v ] = rolling_window_forecast( mdl,insmpl,outofsmpl,0 );
%We then  do rolling-window forecasts using different r.h.s. variables. 
%                                                     right hand side
[E,V,logL]=infer(estmdl,retOtC);
%We store all the variables in a matrix below
RVmat = [rv , V]; 

%Choose the length of the insample period
% insmpl=504;  2013 + 504 = 2517
inLen = length(dcIn);
%Initialize a metrix to collect the forecasts:
res=[];
%Write a loop to compute the forecasts:
for i=1:size(RVmat,2)% i = 1,2,...9
   res=[res HAR_frcst( rv,RVmat(:,i),inLen )];
end
%Note that for the matrix res, all the odd columns store the true value
%RV, and all the even colums store the forecasts from different models

% h) For the out-of-sample period use the ex-post daily annualized RV as a
% proxy of the forcasting target (true volatility) and compare the
% forecasting performances of the ARMA-GARCH and the HAR-RV models. 
% You should compute standard forecasting statistics, conduct appropriate
% tests and provide conclusions.

%Now to evaluate our forecasting performance for each estimator, we can
%construct MSPE, MAPE and QLIKE measures on the big res matrix:
HAR_result=HAR_eval(res);
%Store the results in a nice table
HARtable=array2table(HAR_result, 'RowNames',RVnames,'VariableNames',...
    {'MAPE','MSPE','QLIKE'});
disp(HARtable);

   %                   MAPE          MSPE        QLIKE 
   %                _________    __________    _______

    % RV             0.0086515    0.00015263    0.30441
    % GJR-GARCH-t     0.011787    0.00027487     0.5434
    
%Now, we can use a modified DM test to test if one model significantly
%outperforms another. We will be comparing HAR-RV against GJR-GARCH(1,1)
% To perform the modified DM test, we firstly
%compute the forecasting errors:
evec=res(:,2:2:end)-res(:,1:2:end);
%evec contains the forecasting error computed by forecast-true for each
%series. Note that the 1st series is the RV benchmark. 

%Initialize a matrix to store all results
DMresult = zeros(1,2);
 [DM,p_value] = dmtest_modified(evec(:,1), evec(:,2));
 DMresult(1,:)=[DM p_value];
% end
DMtbl=array2table(DMresult,'RowNames',{'GARCH vs RV'},'VariableNames',{'TestStat','p_value'});
disp(DMtbl);
%From the table we can see that all TestStats are negative, which suggests
%that the HAR-RK model outperform HAR models with other r.h.s. variables.
%However this difference in terms of MSPE is not significantly different
%from zero since the p-values are all very high, except from HAR-GARCH as
%we see a clear rejection. This indicates that the HAR-GARCH model performs
%significantly worse than HAR-RK in predicting RK with a much larger
%forecasting error.

 %                  TestStat     p_value  
 %                  ________    __________

  %  GARCH vs RV    -6.0161     3.4415e-09


 