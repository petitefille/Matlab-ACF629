%Set up and load the dataset
%clear the pre-existing variables
clearvars;
%set working directory (Use your own path instead! Or browse to where the
%data is stored)
%cd('C:\Users\Eier\Documents\Lancaster\ACF609\Lab4')
%load an excel file into MATLAB as table
dmat=readtable('StockIndexFXDATA.xlsx');
data=dlog(dmat{:,13}); % log returns
T=length(data); % 7224


%Picking the best model from ARMA(0,0) to ARMA(pmax,qmax) based on information
%criteria

%initiallizing information criteria mat
pmax=1;
qmax=1;
AICmat=zeros(pmax+1,qmax+1); % why 2x2 matrix ??=> because starting from pmax
                             % = 0 and qmax = 0
BICmat=zeros(pmax+1,qmax+1);
HQICmat=zeros(pmax+1,qmax+1);

%Compute the log-likelihood and information criteria for the models considered
%Note that I started with ARMA(0,0) and loop over all the models to
%ARMA(pmax,qmax).
for i=1:pmax+1 % pmax = 0,1
    for j=1:qmax+1 % qmax = 0,1
        Mdltemp=arima(i-1,0,j-1);
        [~,~,logL,~] = estimate(Mdltemp,data,'Display','off');
             
        %You can also use the aicbic function to compute AIC or BIC.
        % The number of parameters in a model is p + q + 1 (for the AR and MA coefficients, and constant term)
        AICmat(i,j)=-2*logL/T+2*(-1+i+j)/T; 
        % i+j-1 
        BICmat(i,j)=-2*logL/T+(-1+i+j)*log(T)/T;
        HQICmat(i,j)=-2*logL/T+2*(-1+i+j)*log(log(T))/T;
    end
end
%Find the location of the minimum information crtieria. Pay close attention
%to how the find() function is used in here.
[paic,qaic]=find(AICmat==min(min(AICmat)));  % her returneres indeksen
                                            % til minste element i AICmat 
                                            
% A er en matrise: min(A) returnerer en vektor a med minste element av hver kolonne 
% a er en vektor: min(a) returnerer minste element av vektoren
[pbic,qbic]=find(BICmat==min(min(BICmat)));
[phqic,qhqic]=find(HQICmat==min(min(HQICmat)));

%Choose best model using BIC for example
pbest=pbic-1; % 1
qbest=qbic-1; % 1

[P,Q] = ARMA_optimal(data,pmax,qmax,'BIC');% P = 1, Q = 1
%Model Evaluation
%Estimate the best model and compute the residual
Mdlbest=arima(pbest,0,qbest);
[EstMdlbest,EstParamCov,logL1,info] = estimate(Mdlbest,data);
%EstMdlbest: model containing parameter estimates

    % Description: "ARIMA(1,0,1) Model (Gaussian Distribution)"
    % Distribution: Name = "Gaussian"
               % P: 1
               % D: 0
               % Q: 1
        % Constant: -1.51587e-05
        %      AR: {0.781633} at lag [1]
        %     SAR: {}
        %      MA: {-0.879475} at lag [1]
        %     SMA: {}
    % Seasonality: 0
    %       Beta: [1×0]
    %    Variance: 0.00383089

%EstParamCov: Variance-covariance matrix of maximum likelihood estimates of model parameters known to the 
% optimizer, returned as a matrix.
% logL1: Optimized loglikelihood objective function value, returned as a scalar
% info: summary information
%ARIMA(1,0,1) Model (Gaussian Distribution):
 
%                   Value       StandardError    TStatistic    PValue 
%                ___________    _____________    __________    _______

%    Constant    -1.5159e-05     0.00010082       -0.15035     0.88049
%    AR{1}           0.78163       0.017726         44.095           0
%    MA{1}          -0.87948       0.013985        -62.887           0
%    Variance      0.0038309     3.6321e-05         105.47           0

[E,V,logL2]=infer(EstMdlbest,data);
% E: residuals
% V: conditional variances
% logL: the loglikelihood objective function values
% ARIMA(1,0,1) Model (Gaussian Distribution):
 
%                   Value       StandardError    TStatistic    PValue 
%                ___________    _____________    __________    _______

%    Constant    -1.5159e-05     0.00010082       -0.15035     0.88049
%    AR{1}           0.78163       0.017726         44.095           0
%    MA{1}          -0.87948       0.013985        -62.887           0
%    Variance      0.0038309     3.6321e-05         105.47           0

%Let's look at some graphs first...
fig = figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1);
plot(E);
title('Line plot of the residual');
subplot(2,2,2);
histogram(E);
title('Histogram of the residual');
subplot(2,2,3);
%quantile-quantile plot against a standardized residual. Therefore you need
%to standardize the residual first using the zscore() function
qqplot(zscore(E));
title('Quantile-Quantile plot of the standardized residual');
subplot(2,2,4);
%This is exactly the correlogram codes that I used in Computer Lab 1...
hold on;
lag=100;
[acf1,~,bound]=autocorr(E,lag); % 101 x 1 (from la 0 to lag 100)
[acf2,~,bound1]=autocorr(E.^2,lag); % bound and bound1 are the same
p1=stem(1:1:lag,acf1(2:lag+1),'filled'); % ploting acf in blue and filled
p2=stem(1:1:lag,acf2(2:lag+1)); %ploting acf^2 not filled in red
p=plot(1:1:lag,[ones(lag,1).*bound1(1) ones(lag,1).*bound1(2)]); %ploting confidence bounds
hold off;
%Graph settings
p1.Color='r';
p1.MarkerSize=4; % acf1: red and filled
p2.Color='b';
p2.MarkerSize=4; % acf2: blue and not filled
p(1).LineWidth=2; % upper confidence bound
p(1).LineStyle=':'; 
p(1).Color='k'; % in black 
p(2).LineWidth=2;
p(2).LineStyle=':';
p(2).Color='k'; % lower confidence bound
title('Correlogram of the residual and squared residual');
legend('Residual', 'Squared residual','95% Confidence Bounds', 'Location','southoutside','Orientation','horizontal');

%Now..test for normality, we use the Jacque-Bera test. Prepare a nice table
%to present your  results!
%Set the significance level of the test
alpha=[0.1,0.05,0.01];
[h1,p1, JBstat1,critval1]=jbtest(E,alpha(1)); % non-normal
[h2,p2, JBstat2,critval2]=jbtest(E,alpha(2)); % non-normal
[h3,p3, JBstat3,critval3]=jbtest(E,alpha(3)); % non-normal
jbtable=[JBstat1 critval1 critval2  critval3 p1 ];
jbout=array2table(jbtable,'VariableNames',  {'JBstat','CV10', 'CV5', 'CV1','p_value' });
disp('Jacque-Bera test on the residual');
disp(jbout);
% Jacque-Bera test on the residual
%    JBstat     CV10      CV5       CV1      p_value
%    ______    ______    ______    ______    _______

%    9148.3    4.5826    5.9843    9.3521     0.001 


%Ljung-Box test on the residual (for up to 50 Lags) to test for
%autocorrelation
Lags=50;
%Note that we need to adjust the degree of freedom in the Ljung-Box test on the
%residuals. Read help lbtest for details.
[~,pValue,stat,~] = lbqtest(E,'dof',-pbest-qbest+3:Lags-pbest-qbest,'Lags',3:Lags);
% nr
% Expected a string scalar or character vector for the parameter name.
 lbqtable=[ (pbest+qbest+1:Lags)' stat'  pValue' ];
lbqout=array2table(lbqtable,'VariableNames',  {'Lags','LBstat','p_value' });
disp('Ljung-Box test on the residual');
disp(lbqout);

%ARCH LM test to check for potential heteroscedasticity in the  error
%term. We only check for the first lag.
  [h,pValue,stat,cValue] = archtest(E,'Lags',1);
  % Expected a string scalar or character vector for the parameter name.
%You should be able to check the test results and produce a table...