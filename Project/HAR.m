clearvars; clc;
%Load the dataset of RV measures
cd('C:\Users\Eier\Documents\Lancaster\ACF609\Lab7')
load('RVtbl.mat');
% RVtble: 2517 x 10 
%Plot the RV measures. 
%plot(RVtbl{:,2:10}); Fig 1 

%Get the names of the series
RVnames=RVtbl.Properties.VariableNames(2:10);
%We choose the benchmark model to be RK as the l.h.s. of a HAR model
%                                              left hand side 
% {'tickRV'}    {'RV'}    {'SubRV'}    {'TSRV'}    {'RK'}    {'preRV'}    {'RBV'}    {'preRBV'}    {'GARCH'}
RK=RVtbl.RK; % 2517 x 1 

%We then  do rolling-window forecasts using different r.h.s. variables. 
%                                                     right hand side 
%We store all the variables in a matrix below
RVmat=RVtbl{:,2:10}; % 2517 x 9


%Choose the length of the insample period
insmpl=504; % 2013 + 504 = 2517
%Initialize a metrix to collect the forecasts:
res=[];
%Write a loop to compute the forecasts:
for i=1:size(RVmat,2)% i = 1,2,...9
   res=[res HAR_frcst( RK,RVmat(:,i),insmpl )]; % size: 2013 x 18 
   %     res er hva vi har lagret før i res
end
% res 
%Note that for the matrix res, all the odd columns store the true value
%RK, and all the even colums store the forecasts from different models, in
%the same order as stored in the RVtbl.mat

%Now to evaluate our forecasting performance for each estimator, we can
%construct MSPE, MAPE and QLIKE measures on the big res matrix:
HAR_result=HAR_eval(res);% 9 x 3 matrise
%Store the results in a nice table
HARtable=array2table(HAR_result, 'RowNames',RVnames,'VariableNames',...
    {'MAPE','MSPE','QLIKE'});
disp(HARtable);

 %                 MAPE          MSPE        QLIKE 
              __________    __________    _______

 %   tickRV    8.0917e-05    4.7462e-08    0.28455
 %   RV         8.029e-05    4.8377e-08    0.28304
 %   SubRV     8.0846e-05    4.7185e-08    0.28137
 %   TSRV      8.1885e-05    5.1172e-08     0.2944
 %   RK        7.9252e-05    4.4487e-08    0.27677 % less noisy benchmark
 %    preRV     7.9483e-05    4.4965e-08    0.27553 %prerv qlike lowest 
 %   RBV       8.1834e-05    5.1018e-08    0.29021
 %   preRBV    7.9601e-05    4.5145e-08    0.28097
 %   GARCH     9.2406e-05    5.9163e-08    0.34501 % more dependence
 %   between days 
 % IC : 2 IC says RK is best while preRV is best (IC:1)
%From the table we see that different loss functions lead to differnt
%conclusions for the best model, so according to MAPE and MSPE, RK performs the
%best, but QLIKE suggests that preRV performs the best. Also, all three
%loss functions suggest that GARCH measure performs the worst, which is not
%surprising because GARCH measure do not take intraday price information
%into account

%Now, we can use a modified DM test to test if one model significantly
%outperforms another. We will be comparing HAR-RK against all possible
%alternatives in the dataset. To perform the modified DM test, we firstly
%compute the forecasting errors:
evec=res(:,2:2:end)-res(:,1:2:end);
     % res[ kol2 kol4 ...kol8] - res[kol1 kol3 .. kol9] = forecasted - true
% size(evec) 2013 x 9   
% evec[5] = RK = true != 0
%evec contains the forecasting error computed by forecast-true for each
%series. Note that the 5th series is the RK benchmark. To make our
%programme simpler, we rearrange the order of the collums and move RK
%the first column of evec:
evec=evec(:,[5 1:4 6:9]); % flytter RK (true) til første kolonne

%Initialize a matrix to store all results
DMresult=zeros(8,2); % compare 8 RV values to RK (true)
for i=1:8
    %We are testing HAR-RK against all other alternatives. Therefore the
    %first entries will always be the loss from HAR-RK
    % e1: RK 
    % e2: andre 8 RV modeller
 [DM,p_value] = dmtest_modified(evec(:,1), evec(:,i+1));
 DMresult(i,:)=[DM p_value];
end
DMtbl=array2table(DMresult,'RowNames',{'tickRV vs RK'...
    'RV vs RK' 'SubRV vs RK' 'TSRV vs RK' 'preRV vs RK' 'RBV vs RK'...
    'preRBV vs RK' 'GARCH vs RK'},'VariableNames',{'TestStat','p_value'});
disp(DMtbl);
%From the table we can see that all TestStats are negative, which suggests
%that the HAR-RK model outperform HAR models with other r.h.s. variables.
%However this difference in terms of MSPE is not significantly different
%from zero since the p-values are all very high, except from HAR-GARCH as
%we see a clear rejection. This indicates that the HAR-GARCH model performs
%significantly worse than HAR-RK in predicting RK with a much larger
%forecasting error.

%                    TestStat     p_value 
%                    ________    _________

%    tickRV vs RK    -0.81038      0.41782  % not stat sign diff - fail to
                                            % rekect that models are stat
                                            % sign diff
 %   RV vs RK        -0.74277      0.45771
 %   SubRV vs RK     -0.86881      0.38506
%    TSRV vs RK      -0.95268      0.34087
%    preRV vs RK     -0.43457      0.66392       
%    RBV vs RK        -1.0536      0.29217
%    preRBV vs RK    -0.45577       0.6486
%    GARCH vs RK      -2.8285    0.0047226   GARCH is worse significantly
                                             %    diff not designed 
                                             % for tick by tick data
