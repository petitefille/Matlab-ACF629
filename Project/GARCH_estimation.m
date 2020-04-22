clearvars; clc;
%Load the dataset
load('IBM_price.mat');
%Initializing the return vector
% 5034/2 =  2517
ret=zeros(2517,1);
%Get the unique day identifier from the dataset (which is date in this
%case)
% C = unique(A) returns the same data as in A, but with no repetitions. C is in sorted order.
% If A is a table or timetable, then unique returns the unique rows in A in sorted order.
ind=unique(price.date); % 2517 x 1 
%extract price and date vector from the dataset
p=price.price; % 5034 x 1
d=price.date; % 5034 x 1
%Method 1: looping over each day to compute the open-to-close return
for i=1:2517
   temp=p(d==ind(i)); % returns open and close price for each unique day
   % temp = [p1,p2] % returnerer p element når d == ind(i) (2 tilfeller)
    ret(i)=log(temp(end))-log(temp(1));
    %          temp(2) 
    % temp(1) is open price 
    % teemp(2) close price
    % Yt = ln(Pt/Pt-1) = ln(close price/open price) = ln(close price) -
    % ln(open price)
end
%Method 2: use reshape to assign the open and close price into another
%matrix and compute log return.
temp2=reshape(p,2,2517); % 2 x 2517 
% her blir p omgjort til wn 2 x 2517 matrise slik at hver kolonne i temp2
% inneholder open and close price for hver unike dato
% Eks: temp2(:,1) = [open price date 1, close price date 1] = 2 x 1 
ret2=log(temp2(2,:)')-log(temp2(1,:)');
% temp2(2,:)': 2517 x 1 vektor some inneholder alle close price
% temp2(1,:)': 2517 x 1 vektor som inneholder alle open price

%ploting the line graph, histogram (with fit to a normal density),
%autocorrelation for return and autocorrelation for the squared returns.
%Relate this to the stylized facts in the lecture notes!
% Fig 1
subplot(2,2,1);
plot(ret); % log returns 
subplot(2,2,2);
histfit(ret); % histogram/density plot of log returns 
% histfit(data) plots a histogram of values in data using the number of bins equal to the square root of the number of elements in
% data and fits a normal density function.
subplot(2,2,3);
autocorr(ret,1000); % 1000 Lags
subplot(2,2,4);
autocorr(ret.^2,1000); % 1000 Lags - 

%Let us assume that the optimized ARMA structure is ARMA(1,1). Compute
%optimal GARCH lags using a function GARCH_optimal.m. See the function file
%for details.
[gbest,abest,ICmat] = GARCH_optimal(ret, 2, 2,1, 1, 'BIC');
% gbest = 1
% abest = 1
% gmax: 2. Maximum GARCH lags considered in the model 
% amax: 2. Maximum ARCH lags considered in the model
% arpar: 1.Number of AR parameters in the mean equation
% mapar: 1. number of MA parameters in the mean equation
% method: string. Methods to used for choosing the optimal model. Choice:
%'AIC': Akaike Information Criterion
%'BIC': Bayesian Information Criterion
%'HQIC': Hannan-Quinn Information Criterion
%Outputs:
%gbest: optimal GARCH lags = 1
%abest: optimal ARCH lags = 1
%ICmat: Matrix of the chosen information criteria

%Using the optimal GARCH lags 
                                           % garch(1,1) https://se.mathworks.com/help/econ/garch.html
Mdl=arima('ARLags',1','MALags',1,'Variance',garch(gbest,abest));
%Estimate the model
[estmdl, estParamCov,logL, info] =estimate(Mdl,ret);
%Pay attention to the GARCH and ARCH paramters. ARCH+GARCH very close to 1,
%suggesting that the volatility process is very persistent.
 

%    ARIMA(1,0,1) Model (Gaussian Distribution):
 
%                  Value       StandardError    TStatistic      PValue   
%                __________    _____________    __________    ___________

%    Constant    0.00014307     5.6657e-05        2.5252         0.011563
%    AR{1}          0.88452        0.04154        21.293      1.3136e-100
%    MA{1}         -0.84785        0.04756       -17.827       4.3719e-71

 
 
%    GARCH(1,1) Conditional Variance Model (Gaussian Distribution):
% 
%                  Value       StandardError    TStatistic      PValue  
%                __________    _____________    __________    __________

%    Constant    3.0732e-06     7.9474e-07         3.867      0.00011019
 %   GARCH{1}       0.87601       0.012328        71.058               0
 %   ARCH{1}       0.095708       0.009778        9.7881      1.2662e-22

%Inferring the residuals, conditional variance process and the
%log-likelihood from the fitted model
[E,V,logl]=infer(estmdl,ret);
% E: inferred residuals
% V: inferred conditional variances
% logL: loglikelihood objective function values
%Compute fitted standardized error term
z=E./sqrt(V);

%Diagnostic testing
%testing against a standard normal distribution
qqplot(z);
% Jarque Bera test
[h,p,stat,cv]=jbtest(z); % testing standardized residuals
% h: 1 => not normally distributed

%Detecting remaining ARCH effects
[h2,p2,stat2,cv2] = archtest(z,'Lags',1:10);
% h2: 0   0   0   0   0   0   0   0   0   0 0 => no conditional
% heteroscedasticity
% our model captures the heteroscedasticity present in data

% You can also use Ljung-Box tests on z and z.^2 to test for remaining
% autocorrelations...
% does our model capture the dependence of data
% [h3,pValue3,stat3cv3] = lbqtest(E,'dof',-pbest-qbest+3:Lags-pbest-qbest,'Lags',3:Lags);
% [~,pValue,stat,~] = lbqtest(E,'dof',-pbest-qbest+3:Lags-pbest-qbest,'Lags',3:Lags);

%:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%t-distributed GARCH models
%Changing our definition of Mdl slightly
Mdl2=Mdl;% Mdl.Distribution = "Gaussian"
Mdl2.Distribution='t';
%Estimate the model
[estmdl2, EstParamCov2, logL2, info2]=estimate(Mdl2,ret); % uses maximum likelihood to estimate the parameters of the ARIMA(p,D,q) model Mdl given the observed
% univariate time series y. EstMdl is an arima model that stores the results.

 
  %  ARIMA(1,0,1) Model (t Distribution):
 
  %                Value       StandardError    TStatistic      PValue  
%                __________    _____________    __________    __________
%
%    Constant    0.00026862     0.00014264        1.8832        0.059667
%    AR{1}          0.76544        0.12053        6.3507      2.1439e-10
%    MA{1}          -0.7303        0.12887       -5.6671      1.4527e-08
%    DoF             8.3607         1.3452        6.2154      5.1195e-10

 
 
%    GARCH(1,1) Conditional Variance Model (t Distribution):
% 
%                  Value       StandardError    TStatistic      PValue  
%                __________    _____________    __________    __________

 %   Constant    2.9045e-06     9.6994e-07        2.9945       0.0027488
 %   GARCH{1}        0.8727       0.017162        50.851               0
 %   ARCH{1}       0.095745       0.013518        7.0828      1.4127e-12
 %   DoF             8.3607         1.3452        6.2154      5.1195e-10
    
    
[E2,V2,logl2]=infer(estmdl2,ret);
% E2: inferred residuals
% V2: inferred conditional variances
% logL2: loglikelihood objective function values

%Compute fitted standardized error term
z2=E2./sqrt(V2);
%Notice that DoF is around 8, which suggests an excess kurtosis of about
%1.5.
%We need to make sure that our specification using the t-distribution is 
%correct. We can use a Q-Q plot to check it.
%Firstly, we get the estimated DoF and define a t-distribution with the
%estimated DoF
dof=estmdl2.Distribution.DoF; % 8.3607
dist=makedist('tLocationScale',0,1,dof);
% makedist: creates probability distribution name object 
% distribution: tLocationScale
% mu = 0;
% sigma = 1;
% dof = 8.3607

% dist
% t Location-Scale distribution
%        mu =       0
%    sigma =       1
%       nu = 8.36071

%We can then plot the standardized residual accordingly
% to compare the distribution of our standardized residuals and
% the t distribution (0,1,8.35071)
qqplot(z2,dist); % Dif 3
%It is clear that the goodness-of-fit is much better compared to a normal
%assumption.


%Asymmetric GARCH. We consider only EGARCH: https://se.mathworks.com/help/econ/egarch.html
mdl3=arima('ARLags',1,'MALags',1,'Variance',egarch(1,1));
[estmdl3, EstParamCov3,logL3, info3]=estimate(mdl3,ret);  % additionally returns EstParamCov, the variance-covariance matrix 
% associated with estimated parameters, logL, the optimized loglikelihood objective function, and info, a data structure of 
% summary information.

% arima with properties:

%     Description: "ARIMA(1,0,1) Model (Gaussian Distribution)"
%   Distribution: Name = "Gaussian"
%               P: 1
%               D: 0
%               Q: 1
%        Constant: 7.70309e-05
%              AR: {0.912835} at lag [1]
%             SAR: {}
%              MA: {-0.873253} at lag [1]
%             SMA: {}
%     Seasonality: 0
%            Beta: [1×0]
%        Variance: [EGARCH(1,1) Model]
        

[E3,V3,logL3]=infer(estmdl3,ret); 
% E3: inferred residuals
% V3: inferred conditional variances
% logL3: loglikelihood objective function values

%Compute fitted standardized error term
z3=E3./sqrt(V3);
% assessing EGARCH model 

%and GJR-GARCH:
mdl4=arima('ARLags',1,'MALags',1,'Variance',gjr(1,1)); % https://se.mathworks.com/help/econ/gjr-model.html
[estmdl4,EstParamCov4,logL4,info4] =estimate(mdl4,ret); % additionally returns EstParamCov, the variance-covariance matrix associated with estimated 
% parameters, logL, the optimized loglikelihood objective function, and info, a data structure of summary information.

%    ARIMA(1,0,1) Model (Gaussian Distribution):
% 
%                  Value       StandardError    TStatistic      PValue   
%                __________    _____________    __________    ___________
%
%    Constant    9.2223e-05     3.9313e-05        2.3459         0.018983
%    AR{1}          0.90098       0.030602        29.442      1.6162e-190
%    MA{1}          -0.8575       0.037081       -23.125      2.5793e-118

 
 
%    GJR(1,1) Conditional Variance Model (Gaussian Distribution):
 
%                     Value       StandardError    TStatistic      PValue  
%                   __________    _____________    __________    __________
%
%    Constant       3.2189e-06     7.6147e-07        4.2273      2.3656e-05
%    GARCH{1}          0.88137       0.012345        71.392               0
%    ARCH{1}          0.029224       0.011595        2.5204        0.011721
%    Leverage{1}       0.11832       0.017786        6.6527      2.8781e-11


[E4,V4,logl4]=infer(estmdl4,ret);
z4=E4./sqrt(V4); % standardized residuals
% E4: inferred residuals
% V4: inferred conditional variances
% logL4: loglikelihood objective function values
%How to interpret the leverage parameter and how does it relates to the
%volatiltiy asymmetric effect? See lecture notes for details.

%Of course they can also be specified with a t-distirbuted error term.
%Since the GJR-GARCH model nests the GARCH model, we can conduct a
%likelihood-ratio test to see if the extra parameter in the GJR-GARCH model
%improves the likelihood significantly.
 [h3,p3,stat3,cv3] = lratiotest(logl4,logl,1);
 % logl4:  unrestricted model
 % logl: restricted model
 % h3 = 1; p3 = 1.8839e-09 
 % h3 = 1 means unrestricted model (GJR-GARCH) is better
%The test is significant at any conventional significance level, suggesting
%that we should use an asymmetric structure when modelling financial
%returns with GARCH-type structure.

%The estimated volatility series is stored in the V vector. We can plot the
%estimated volatility series from the four different models here:
plot([V V2 V3 V4]);
legend({'GARCH','GARCH-t','EGARCH','GJR-GARCH' },'FontSize',14);
%Note that they hav every similar shapes. However, if you zoom in at the
%financial crisis period (1000) you will see that the asymmetric GARCH models
%produces larger volatility estimates due to the large price drops. 

%Finally, the forecast() function can be used to forecast conditional mean
%and variance simultaneously in the same fashion as forecasting an ARMA
%model. For the general setting, see solutions to the previous workshop.

