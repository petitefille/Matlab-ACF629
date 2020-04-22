
%Set up and load the dataset
%clear the pre-existing variables
clearvars;
%set working directory (Use your own path instead! Or browse to where the
%data is stored)
% cd('E:\Dropbox\Teaching\AcF609-2017\Computer Lab 1');
cd('C:\Users\Eier\Documents\Lancaster\ACF609')
%load an excel file into MATLAB as table
dmat=readtable('StockIndexFXDATA.xlsx');
%extract the numerical part of the data (note that I used CURLY brackets
%to index tables. This returns the ACTUAL elements stored in the table.
%Since the elements are all numerical, the output variable is a numerical
%matrix. If the PARANTHESES are used, such as dmat(:,5:13), then the
%output variable is a TABLE instead of a numerical matrix.
data=dmat{:,5:5}; % type matrise - bruk {} for å konvertere
                   % (:,5:13) = alle rader og kolonne 5-13

%(optional) compute log return
%see dlog.m in your directory for details.
data=dlog(data); % data: 7724x9 matrise 
                 % regne ut log returns for alle de 9 indeksene

%get the size of the data
[rows, cols]=size(data); % 7224 x 9 
%storing the date of the data
datevec=dmat.DATE; % type: date ?
% size(datavec) =  7224 x 1
% datevec(1)= 02.01.1990
% datevec(7224) = 08.09.2017

% storing the name of the series (This case I used the parentheses to get a
% cell array instead of character arrays!)
names=dmat.Properties.VariableNames(5:13);

% names =

 % 1×9 cell array

  % Columns 1 through 7

%     {'DAX30'}    {'DJ30'}    {'NIKKEI225'}    {'FTSE100'}    {'CAC40'}    {'EUROSTOXX50'}    {'SP500'}

%   Columns 8 through 9

 %    {'SMI'}    {'VIX'}




%Computing the statistics
%mean
% M = mean(A) returns the mean of the elements of A along the first array dimension whose size 
% does not equal 1.
% If A is a vector, then mean(A) returns the mean of the elements.
% If A is a matrix, then mean(A) returns a row vector containing the mean
% of each column. 
meandata=mean(data);
% meandata =

   % 1.0e-03 *  % OBS

    % 0.2650    0.2836   -0.0973    0.1535    0.1307    0.1580    0.2662    0.2231   -0.0488

    % mean(data(:,1))= 2.6497e-04
%variance
vardata=var(data);
%maximum and miminum
stat = [meandata; vardata; maxdata;mindata; stddata; skedata;kurdata;qdata];
% vardata =

%    0.0002    0.0001    0.0002    0.0001    0.0002    0.0002    0.0001    0.0001    0.0039

% minimum 
mindata=min(data);
%    -0.0987   -0.0820   -0.1211   -0.0927   -0.0947   -0.0901   -0.0947   -0.0907   -0.3506

% maximum

maxdata = max(data);
%     0.1080    0.1051    0.1323    0.0938    0.1059    0.1044    0.1096    0.1079    0.4960
%standard deviation. Note that std(VIX)^2=var(VIX)
stddata=std(data);
%Computing skewness using the built-in function skewness
%     0.0140    0.0104    0.0148    0.0109    0.0137    0.0133    0.0109    0.0114    0.0627

% skewness = the sample skewness of X
skedata=skewness(data);
%    -0.1613   -0.1664   -0.1352   -0.1258   -0.0633   -0.1192   -0.2493   -0.2464    0.7068


%alternatively, one can write an anonymous function for skewness:
skew=@(x) mean((x-mean(x)).^3)./var(x).^(3/2); % def of skewness calculation
                                               % https://en.wikipedia.org/wiki/Skewness
skedata2=skew(data);
%    -0.1612   -0.1664   -0.1352   -0.1258   -0.0633   -0.1192   -0.2493   -0.2463    0.7067

%or write a proper function 'myskewness'. See the end of the script.
%functions can also be called from a separate function script in the working
%directory, just as the dlog function
skedata3=myskewness(data);
%   -0.1612   -0.1664   -0.1352   -0.1258   -0.0633   -0.1192   -0.2493   -0.2463    0.7067

%kurtosis
% sample kurtosis 
kurdata=kurtosis(data);
% 7.9880   11.7317    8.8594    9.3023    7.7720    8.3973   12.2978    9.8247    7.7093

%quantiles
qvector=[.01 .05 .1 .5 .9 .95 .99]; % vektor 
qdata=quantile(data,qvector);
% -0.0414   -0.0297   -0.0403   -0.0306   -0.0400   -0.0395   -0.0306   -0.0336   -0.1474
%   -0.0224   -0.0160   -0.0237   -0.0166   -0.0218   -0.0212   -0.0169   -0.0174   -0.0924
%   -0.0156   -0.0109   -0.0168   -0.0116   -0.0151   -0.0144   -0.0114   -0.0121   -0.0697
%    0.0004    0.0002         0    0.0000         0    0.0003    0.0002    0.0003   -0.0008
%    0.0153    0.0113    0.0165    0.0116    0.0150    0.0140    0.0114    0.0121    0.0703
%    0.0213    0.0155    0.0225    0.0165    0.0207    0.0202    0.0161    0.0166    0.1019
%    0.0352    0.0281    0.0375    0.0289    0.0350    0.0352    0.0302    0.0295    0.1884

% 60 linjer 

%Prepare a table for the results
%combine the statistics first:
stat=[meandata; vardata; maxdata;mindata; stddata; skedata; kurdata;qdata];
%matrise = [rad1;rad2;...;radn]
% stat =

% meandata    0.0003    0.0003   -0.0001    0.0002    0.0001    0.0002    0.0003    0.0002   -0.0000
% vardata    0.0002    0.0001    0.0002    0.0001    0.0002    0.0002    0.0001    0.0001    0.0039
% maxdata    0.1080    0.1051    0.1323    0.0938    0.1059    0.1044    0.1096    0.1079    0.4960
% mindata   -0.0987   -0.0820   -0.1211   -0.0927   -0.0947   -0.0901   -0.0947   -0.0907   -0.3506
% stddata    0.0140    0.0104    0.0148    0.0109    0.0137    0.0133    0.0109    0.0114    0.0627
% skedata   -0.1613   -0.1664   -0.1352   -0.1258   -0.0633   -0.1192   -0.2493   -0.2464    0.7068
% kurdata    7.9880   11.7317    8.8594    9.3023    7.7720    8.3973   12.2978    9.8247    7.7093
% qdata x_0.01    -0.0414   -0.0297   -0.0403   -0.0306   -0.0400   -0.0395   -0.0306   -0.0336   -0.1474
% x_0.05    -0.0224   -0.0160   -0.0237   -0.0166   -0.0218   -0.0212   -0.0169   -0.0174   -0.0924
% x_0.1   -0.0156   -0.0109   -0.0168   -0.0116   -0.0151   -0.0144   -0.0114   -0.0121   -0.0697
% x_0.5    0.0004    0.0002         0    0.0000         0    0.0003    0.0002    0.0003   -0.0008
% x_0.9    0.0153    0.0113    0.0165    0.0116    0.0150    0.0140    0.0114    0.0121    0.0703
% x_0.95    0.0213    0.0155    0.0225    0.0165    0.0207    0.0202    0.0161    0.0166    0.1019
% x_0.99    0.0352    0.0281    0.0375    0.0289    0.0350    0.0352    0.0302    0.0295    0.1884
    
    
%convert it to table with proper titles
% names er kolonne navn
stattable=array2table(stat,'VariableNames',names,'RowNames', {'Mean' ...
    'Variance' 'Maximum' 'Minimum' 'Std.Dev.' 'Skewness' 'Kurtosis' ...
    'Q(0.01)' 'Q(0.05)' 'Q(0.1)' 'Q(0.5)' 'Q(0.9)' 'Q(0.95)' 'Q(0.99)'});
%display the table to the output window
disp(stattable);
 %                 DAX30          DJ30        NIKKEI225      FTSE100        CAC40       EUROSTOXX50      SP500          SMI            VIX    
 %               __________    __________    ___________    __________    __________    ___________    __________    __________    ___________
%
%    Mean        0.00026497    0.00028358    -9.7259e-05     0.0001535    0.00013072    0.00015797     0.00026623    0.00022313    -4.8778e-05
%    Variance     0.0001971    0.00010879     0.00022007    0.00011842    0.00018661    0.00017775     0.00011988    0.00012897      0.0039257
%    Maximum        0.10797       0.10508        0.13235      0.093843       0.10595       0.10438        0.10957       0.10788        0.49601
%    Minimum      -0.098707     -0.082005       -0.12111     -0.092656     -0.094715     -0.090111      -0.094695     -0.090703       -0.35059
%    Std.Dev.      0.014039       0.01043       0.014835      0.010882      0.013661      0.013332       0.010949      0.011357       0.062655
%    Skewness      -0.16128      -0.16642       -0.13522      -0.12583     -0.063322      -0.11921       -0.24932       -0.2464        0.70682
%    Kurtosis         7.988        11.732         8.8594        9.3023         7.772        8.3973         12.298        9.8247         7.7093
%    Q(0.01)      -0.041434     -0.029672      -0.040325     -0.030615     -0.040019     -0.039456      -0.030607     -0.033551       -0.14741
%    Q(0.05)      -0.022447     -0.016048      -0.023723     -0.016585     -0.021759     -0.021151      -0.016873     -0.017442      -0.092424
%    Q(0.1)        -0.01558      -0.01092      -0.016848     -0.011594     -0.015101     -0.014398      -0.011376     -0.012109      -0.069728
%    Q(0.5)      0.00044861    0.00020085              0    2.6987e-05             0    0.00029752     0.00020698    0.00027285    -0.00083058
%    Q(0.9)         0.01533      0.011273       0.016486      0.011575       0.01497      0.013998       0.011422      0.012142        0.07033
 %   Q(0.95)       0.021327      0.015455       0.022463      0.016506       0.02066      0.020158       0.016085      0.016633        0.10193
%    Q(0.99)       0.035205      0.028124       0.037542       0.02891      0.035024      0.035155       0.030154      0.029478        0.18841


%plot the line graph and histogram of the dataset in multiple plots
%set size of the plots [left bottom width height]
plotsize=[100,100,1200,550]; % radvektor
%create a figure structure named 'fig1'
fig1=figure;
for i=1:9
    %subplot(M,N,i) creates a M-by-N subplot matrix of plots.
    %M: number of rows of plots
    %N: number of columns of plots
    %i: the index of plots to be plotted
    subplot(3,3,i);
    %line plot
    plot(datevec,data(:,i)); % log returns for index i 
    %set titles of the subplots
        title(names(i));  % names of index i   
end
% ------------------------


%set the super title of the whole plot
suptitle('Line graphs of the data');
%set the size of the plot
fig1.Position=plotsize; & ??
%save plot
saveas(fig1, 'Linegraph.jpg');

fig2=figure;
for i=1:9
    subplot(3,3,i);
    %histogram
    histogram(data(:,i));
        title(names(i));
end
suptitle('Histograms of the data');
fig2.Position=plotsize;
saveas(fig2, 'Histogram.jpg');


%Correlogram and partial correlogram up to lag L 
L=36;
fig3=figure;
Q=zeros(L,9); % 36  x 9 matrise
for i=1:9
    %computing autocorrelation and partial autocorrelation
    % [acf,lags,bounds] = autocorr(___) additionally returns the lag numbers that MATLAB® uses 
    % to compute the ACF, and also returns the approximate upper and lower confidence bounds.
    [acf,lag1,bound1]=autocorr(data(:,i),L);
    
    [pacf,lag2,bound2]=parcorr(data(:,i),L);
    %compute Ljung-Box test statistics for later
    Q(:,i)=length(data)*(length(data)+2)*cumsum((acf(2:L+1).^2)./(length(data)-lag1(2:L+1)));
    
    subplot(3,3,i);    
    hold on; %adding new plots to the existing plot   
    p1=stem(1:1:L,acf(2:L+1),'filled'); %ploting acf
    p2=stem(1:1:L,pacf(2:L+1)); %ploting pacf
    p=plot(1:1:L,[ones(L,1).*bound1(1) ones(L,1).*bound1(2)]); %ploting confidence bounds
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
    title(names(i)); 
     
        %legend settings
    if i==8
        legend('ACF', 'PACF','95% Confidence Bounds', 'Location','southoutside','Orientation','horizontal');
    end
end
suptitle('Correlogram and partial correlogram of the data');
fig3.Position=plotsize;
saveas(fig3, 'Corgram.jpg');


% ----------------------------------------------------------
%Jacque-Bera test for normality
% H0: Normal
% https://en.wikipedia.org/wiki/Jarque–Bera_test
JB=(length(data))/6*(skedata.^2+0.25*(kurdata-3).^2); % data: log returns
 % 1.0e+04 *

  %   0.7520    2.2982    1.0356    1.1975    0.6859    0.8785    2.6096    1.4093    0.7277
  
%  JB(1)

%  7.5202e+03
% 2 degrees of freedom

JBpvalue=chi2cdf(JB,2,'upper'); % p-verdi
%  0     0     0     0     0     0     0     0     0
%Ljung-Box test for correlation (pvalue computed using the Q obtained from the
%correlogram)
% https://en.wikipedia.org/wiki/Jarque–Bera_test
% P = chi2cdf(X,V,'upper') returns the upper tail probability of 
 %    the chi-square distribution.
 % size(Q) 36 x 9
% https://en.wikipedia.org/wiki/Ljung–Box_test 
% B = repmat(A,r1,...,rN) specifies a list of scalars, r1,..,rN, that describes how copies of 
% A are arranged in each dimension. When A has N dimensions, the size of 
% B is size(A).*[r1...rN]. For example, repmat([1 2; 3 4],2,3) returns a 4-by-6 matrix.
% A'= transpose(A)
LBpvalue=chi2cdf(Q,repmat((1:1:L)',1,9),'upper'); % 1 36x1 kollone repetert 9 ganger 

% repmat((1:1:L)',1,9)

% ans =

%      1     1     1     1     1     1     1     1     1
%      2     2     2     2     2     2     2     2     2
%      3     3     3     3     3     3     3     3     3
% ...
%    34    34    34    34    34    34    34    34    34
%    35    35    35    35    35    35    35    35    35
%    36    36    36    36    36    36    36    36    36

% LBpvalue =

%    0.5771    0.0000    0.0276    0.3524    0.4996    0.2971    0.0000    0.0033    0.0000
 %   0.1503    0.0000    0.0005    0.0012    0.0298    0.0111    0.0000    0.0001    0.0000
 
   % 0.0197    0.0000    0.0077    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000
  %  0.0126    0.0000    0.0093    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000
% prepare a table to present results
teststat=[ JB;Q]; % ikke P verdi
% 37 x 9 
% JB tester for normalitet: 1 x 9 
% Q : LB verdi (correlation up to lag l) : 36 x 9 
pstat=[ JBpvalue; LBpvalue]; % 37 x 9 
testtable=cell(L+1,10); % 37 x 10 
% A cell array is a data type with indexed data containers called cells, where each cell 
% can contain any type of data
%Use * to indicate significance: ***: significant at 1%, **: significant at
%5%, *: significant at 10%

%Try to figure out what I did here!
for i=1:L+1 % 1 - 37 (antall rader) => for hver rad (1)
    if i==1 % første rad 
        testtable(i,1)={'JB stat'}; % cell[1,1] navn 
    else
        testtable(i,1)={strcat('LB(', num2str(i-1,'%5.0f'),')')};% navn på cellen 
            % strcat: concanates strings horizontally 
    end
    for j=2:10  % for hver kolonne (2) 
        %  p value <= alpha: reject H0
        if pstat(i,j-1)<=0.01  % H0 blir forkastet hvis alpha = 0.01 
            testtable(i,j)={strcat(num2str(teststat(i,j-1),'%5.2f'),'***')};% test verdi og ikke p value
                                                                            % settes inn 
        elseif pstat(i,j-1)<=0.05
              testtable(i,j)={strcat(num2str(teststat(i,j-1),'%5.2f'),'**')};
        elseif pstat(i,j)<=0.1
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
% cell2table: converts cell array to table
disp(testtable); % Display value of variable  (disp(X))

% skew=@(x) mean((x-mean(x)).^3)./var(x).^(3/2);
%Functions attached at the end of the script
function [ret]=myskewness(x)
ret=mean((x-mean(x)).^3)./var(x).^(3/2);
end

function [ret]=myskewness(x)
ret=mean((x-mean(x)).^3)./var(x).^(3/2);
end

% skew=@(x) mean((x-mean(x)).^3)./var(x).^(3/2); % def of skewness calculation
                                               % https://en.wikipedia.org/wiki/Skewness
                                               
% where there is a star: H0 is rejected !  -> correlated: meaning?                                              
%no star: not correlated
%                DAX30            DJ30           NIKKEI225         FTSE100          CAC40        EUROSTOXX50         SP500             SMI             VIX     
%              ____________    _____________    _____________    _____________    ____________    ____________    _____________    _____________    ____________
%
%   JB stat    '7520.19***'    '22982.15***'    '10356.09***'    '11974.59***'    '6859.22***'    '8785.46***'    '26096.17***'    '14092.76***'    '7277.04***'
%   LB(1)      '0.31*'         '16.72***'       '4.85**'         '0.86'           '0.46'          '1.09*'         '23.56***'       '8.64***'        '45.52***'  
%   LB(2)      '3.79*'         '25.22***'       '15.16***'       '13.49***'       '7.02**'        '9.01**'        '31.19***'       '19.69***'       '73.90***'  
%   LB(3)      '6.52*'         '25.54***'       '16.68***'       '37.44***'       '25.51***'      '25.57***'      '31.62***'       '26.65***'       '88.46***'  
%   LB(4)      '8.57*'         '25.89***'       '16.73***'       '48.99***'       '27.74***'      '29.11***'      '32.14***'       '35.68***'       '97.46***'  
%   LB(5)      '15.35***'      '29.14***'       '16.74***'       '64.48***'       '45.89***'      '46.42***'      '37.53***'       '59.21***'       '100.10***' 
%   LB(6)      '18.23***'      '30.12***'       '19.32***'       '70.51***'       '46.33***'      '47.17***'      '38.50***'       '65.43***'       '107.77***' 
%   LB(7)      '18.30**'       '37.00***'       '19.41***'       '70.57***'       '46.65***'      '47.39***'      '47.63***'       '68.61***'       '116.69***' 
%   LB(8)      '19.33**'       '41.27***'       '19.52**'        '74.79***'       '48.74***'      '52.18***'      '52.84***'       '72.83***'       '117.46***' 
%   LB(9)      '19.48**'       '41.95***'       '20.55**'        '74.98***'       '49.97***'      '53.31***'      '54.38***'       '72.94***'       '118.71***' 
%   LB(10)     '19.52**'       '46.07***'       '22.23**'        '75.31***'       '50.46***'      '53.35***'      '59.02***'       '73.51***'       '137.47***' 
%   LB(11)     '22.16**'       '46.07***'       '22.42**'        '75.44***'       '50.61***'      '53.66***'      '59.04***'       '74.79***'       '142.86***' 
%   LB(12)     '22.20**'       '46.19***'       '22.63**'        '75.67***'       '50.67***'      '53.66***'      '59.06***'       '74.89***'       '143.40***' 
%   LB(13)     '22.54**'       '52.35***'       '23.47**'        '76.13***'       '50.81***'      '53.97***'      '66.60***'       '78.17***'       '143.50***' 
%   LB(14)     '24.39**'       '53.92***'       '23.67*'         '77.58***'       '50.95***'      '54.01***'      '67.43***'       '78.21***'       '146.08***' 
%   LB(15)     '24.43*'        '55.64***'       '23.72*'         '77.71***'       '51.88***'      '54.18***'      '70.62***'       '78.91***'       '148.06***' 
%   LB(16)     '26.03*'        '56.33***'       '26.47**'        '81.24***'       '57.28***'      '60.78***'      '72.46***'       '78.93***'       '152.90***' 
%   LB(17)     '26.03*'        '57.97***'       '26.55*'         '81.40***'       '57.44***'      '60.92***'      '75.47***'       '79.48***'       '152.91***' 
%   LB(18)     '30.23**'       '69.57***'       '26.93*'         '93.24***'       '61.13***'      '65.60***'      '85.31***'       '81.70***'       '153.27***' 
%   LB(19)     '31.04**'       '69.71***'       '26.93*'         '95.70***'       '62.37***'      '68.13***'      '85.31***'       '81.71***'       '154.36***' 
%   LB(20)     '31.61**'       '69.96***'       '28.59*'         '95.72***'       '63.06***'      '69.28***'      '85.32***'       '86.10***'       '156.68***' 
%   LB(21)     '31.87*'        '69.96***'       '33.05**'        '98.36***'       '63.06***'      '69.54***'      '85.32***'       '86.73***'       '157.32***' 
%   LB(22)     '31.93*'        '70.89***'       '36.80**'        '98.42***'       '63.18***'      '69.59***'      '85.85***'       '86.84***'       '158.06***' 
%   LB(23)     '32.05*'        '72.36***'       '36.80**'        '100.80***'      '63.99***'      '70.21***'      '87.06***'       '86.85***'       '158.79***' 
%   LB(24)     '32.67*'        '75.45***'       '38.38**'        '101.67***'      '64.74***'      '71.26***'      '90.38***'       '86.87***'       '160.64***' 
%   LB(25)     '34.55*'        '75.63***'       '39.43**'        '109.96***'      '65.74***'      '72.30***'      '90.82***'       '87.02***'       '165.05***' 
%   LB(26)     '34.87*'        '81.12***'       '39.47**'        '114.37***'      '66.01***'      '72.35***'      '92.19***'       '87.47***'       '166.42***' 
%   LB(27)     '37.90*'        '81.91***'       '44.04**'        '118.28***'      '69.11***'      '76.14***'      '92.73***'       '87.99***'       '168.82***' 
%   LB(28)     '37.95*'        '84.44***'       '44.33**'        '118.29***'      '69.11***'      '76.48***'      '97.46***'       '87.99***'       '169.99***' 
%   LB(29)     '42.81**'       '84.78***'       '48.38**'        '120.52***'      '73.41***'      '85.56***'      '97.60***'       '89.31***'       '174.78***' 
%   LB(30)     '44.68**'       '85.34***'       '49.69**'        '121.10***'      '78.06***'      '87.62***'      '98.23***'       '101.87***'      '174.82***' 
%   LB(31)     '44.69*'        '86.05***'       '50.09**'        '133.58***'      '79.95***'      '90.04***'      '98.47***'       '103.98***'      '175.02***' 
%   LB(32)     '46.05*'        '86.13***'       '54.85***'       '134.99***'      '82.61***'      '96.35***'      '98.76***'       '104.72***'      '175.02***' 
%   LB(33)     '46.30*'        '86.13***'       '57.47***'       '137.87***'      '83.69***'      '97.76***'      '100.01***'      '104.72***'      '175.08***' 
%   LB(34)     '52.67**'       '89.45***'       '58.36***'       '142.71***'      '87.70***'      '104.44***'     '107.06***'      '105.48***'      '182.50***' 
%   LB(35)     '54.32**'       '102.63***'      '58.48***'       '145.26***'      '89.72***'      '105.81***'     '121.04***'      '106.44***'      '182.97***' 
%   LB(36)     '57.59**'       '103.01***'      '58.91***'       '147.79***'      '90.87***'      '108.23***'     '121.04***'      '106.87***'      '186.33***' 

                                               
                                               