clearvars; clc;
%In  this script we simulate an ARMA model first and attempt to estimate
%the pre-defined parameters. We use the GARMA.m function that we define
%previously.
c=0.3;
phi=0.7;
theta=0.2;
sigmasq=2;
tvec=[c phi theta sigmasq];
y = GARMA( 10000, 1, tvec, 1); % ??


%Method 1: we use the built-in arima structure to estimate the parameters.
%The outputs are directly printed to the command window. 
Mdl=arima(1,0,1); % ARMA(1,1) process
[EstMdl,EstParamCov,logL,info] = estimate(Mdl,y);% y må være en vektor

ARIMA(1,0,1) Model (Gaussian Distribution):
 
                 Value     StandardError    TStatistic      PValue  
                _______    _____________    __________    __________

    Constant    0.27901       0.019286        14.467      1.9557e-47
    AR{1}       0.70874      0.0089552        79.143               0
    MA{1}       0.19411       0.012546        15.472      5.3493e-54
    Variance     2.0338       0.028219         72.07               0
%The parameter estimates are stored in the structure EstMdl. To retrieve
%them, use the following command:
bhat=[EstMdl.Constant; EstMdl.AR{1};  EstMdl.MA{1}; EstMdl.Variance   ];

bhat =

    0.2790 % c
    0.7087 % phi
    0.1941 % theta 
    2.0338 % var of epsilon ??
%The standard errors can be retrieved from the EstParamCov matrix

se=sqrt(diag(EstParamCov)); % EstParamCov is Cov(parameters)
    0.0193 SE(c)
    0.0090 SE(phi)
    0.0125 SE(theta)
    0.0282 SE(var of epsilon)
EstParamCov =

   1.0e-03 *

    0.3720   -0.0809    0.0723   -0.0159
   -0.0809    0.0802   -0.0684    0.0022
    0.0723   -0.0684    0.1574   -0.0002
   -0.0159    0.0022   -0.0002    0.7963

%The log-likelihood function is stored in the variable logL
%Note that we can construct t-statistics and p-values ourselves:
tstat=bhat./se; % tstat = estimated value of parameters/ SE(parameters)
% are you standardizing ?
% distribution
14.4671
79.1427
15.4721
72.0702
pvalue=(1-normcdf(abs(tstat)))*2; % all 0 : what is H0? => H0 rejected
%Prepare a table for the estimation results ourselves. Observe that the
%estimated parmaeters are very close to their corresponding true values.
armadata=[bhat se tstat pvalue ];
armaout=array2table(armadata, 'RowNames', {'c' 'phi' 'theta' 'sigmasq'}, 'VariableNames', {'bhat' 'SE' 'tstat' 'pvalue'});
disp(armaout);

                bhat         SE        tstat     pvalue
               _______    _________    ______    ______

    c          0.27901     0.019286    14.467      0   
    phi        0.70874    0.0089552    79.143      0   
    theta      0.19411     0.012546    15.472      0   
    sigmasq     2.0338     0.028219     72.07      0  

%Method 2: manually programme the likelihood function and optimize it!
%You need the log-likelihood function defined in LL_ARMA11.m

%Set some initial parameters
x0=[0 0 0 var(y) ]; % LL_ARMA11(c,phi,theta,sigmasq,y)% 0 0 0 1
%Define an anonymous function that returns the value of NEGATIVE log-likelihood when
%the input is the parameter vector. 
ll=@(x) -LL_ARMA11(x(1), x(2), x(3), x(4), y);
%We use the fminunc function to minimize the negative log-likelihood based
%on some initial guesses of parameters x0. The fminunc will attempt to
%change the value of the parameter vector until the negative log-likelihood is
%minimized, or equivalently, the log-likelihood is maximized.

%Outputs of fminunc: x: the optimized parameters, fval: the function value,
%exitflag: the exit status of the optimizer, output: the detailed output of
%the optimizer. (grad: the estimated gradient evaluated at the estimated
%parameter value. hessian: the estimated Hessian matrix.)
[x,fval,exitflag,output,grad,hessian]=fminunc(ll,x0);
%We can also generate a table to present our estimated parameters along
%with their significance
%computing the standard errors using the inverse of the Hessian matrix
se2=sqrt(diag(inv(hessian)))';
%conducting simple tests for parameters
tstat2=x./se2;
pvalue2=(1-normcdf(abs(tstat2)))*2;
%Prepare a table for the estimation results
armadata2=[x' se2' tstat2' pvalue2' ];
armaout2=array2table(armadata2, 'RowNames', {'c' 'phi' 'theta' 'sigmasq'}, 'VariableNames', {'bhat' 'SE' 'tstat' 'pvalue'});
disp(armaout2);
%The (minimized) negative log-likelihood function is stored in the variable
%fval. Also notice that the estimated parameters from Method 2 is very close to
%those estimated by Method 1. However, the standard errors estimated from Method 2 is
%slightly different from Method 1 due to different variance-covariance matrix
%estimation procedures used in the two methods!