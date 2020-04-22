clearvars;  % clears variables from memory
clc; % Remove items from workspace, freeing up system memory
cd('C:\Users\Eier\Documents\Lancaster\ACF609\Lab2\Computer Lab 2 Solution')
%Method 1: simulation manually
%Step 1: set some values for the parameters and the number of observations
T=10000; % antall observasjoner (tidspunkter)
c=0.3;
phi=0.7;
theta=0.4;
sigmasq=0.1; % Var(epsilon)

epsilon=randn(T,1)*sqrt(sigmasq);
% randn(T,1) - standard normally distributed random numbers
% now epsilon are normally distributed random variables all with mean 0 and
% var = 0.1 (and sd = srqt(0.1))

y=zeros(T,1); % y: 10 000 x 1 (vector)

y(1)=c/(1-phi); % set initial value of y (= y_1) equal to its mean/ expected value  
% matlab: first index starts at 1 (and not at 0)
for i=2:T % = 10 0000
   y(i)=c+phi*y(i-1)+theta*epsilon(i-1)+epsilon(i);     
end

%Method 2: simulate via the arima structure in MatLab

model=arima('Constant',c,'AR',phi,'MA',theta,'Variance',sigmasq);
y2=simulate(model,10000);


%Testing that the simulation works by checking the unconditional mean,
%variance and autocorrelation:

% disp(): Display value of variable
disp('Theoretical unconditional mean: E(y)=c/(1-phi)=1.0000'); 
% Theoretical unconditional mean: E(y) = c/(1-phi) = 1.000
disp(mean([y y2])); % y, y2: 10 000 x 1 => [y y2] = 10000 x 2 matrise
%  1.0022    0.9969  - y is more accurate than y

% mean() If A is a matrix, then mean(A) returns a row vector containing the mean of each column.

% test = ((1 + 2*phi*theta + theta^2)*sigmasq)/(1-phi^2)= 0.3373
disp('Theoretical unconditional variance: Var(y)=sigmasq*(1+(phi+theta)^2/(1-phi^2))=0.4333');%??
% Theoretical unconditional variance: Var(y)=sigmas1*(1 + (phi + theta)^2/(1-phi^2))=0.4333
% sigmasq*(1+(phi+theta)^2/(1-phi^2)) = 0.3373
disp(var([y y2]));
% 0.3459    0.3363  - y2 is more accurate than y

%Compute the theoretical autocorrelation function for the simulated ARMA model, and
%compare it against the realized values...

rho1=(phi+theta)*(1+phi*theta)/(1+2*phi*theta+theta^2); % theoretical value= 0.8186
% != rho0 = 1
rhovec=(rho1*phi.^[0:9])'; % size: 10 x 1 : gamma1 - gamma10

% Realized values: 
% Method 1:
acfy=autocorr(y,10); % 10 is nr of lags (size: 11 x 1) (includes rho0 = 1)
% Corr(Y_t, Y_t-i) for i = 0 - 10
% Method 2: 
acfy2=autocorr(y2,10); % 11 x 1 matrise 
%Generate a plot to compare the autocorrelations. Clearly the theoretical
%values are very close to the realized values.
% theoretical values vs realized values (simulated by me) meaning + which are which 
plot([rhovec acfy(2:end) acfy2(2:end)]); % Fig 1
legend({'Theoretical autocorrelation' 'Autocorrelation of ARMA(1,1) simulated from Method 1'...
'Autocorrelation of ARMA(1,1) simulated from Method 2' }); 

%Write a function that simulates an M-by-N matrix of independent ARMA(1,1) models 
%with each column an independent ARMA(1,1) model. See function file GARMA.m
%Define the output tvec
tvec=[c phi theta sigmasq]; % vector: 1 x 4 
%We try both options to simulate the ARMA process, and print the time used
%to the command window. It is clear that option 1 is much faster ( method
% 1 is much faster than method 2)
%since we program in vectors and compute N processes simultaneously.
% N: nr of columns = nr of different ARMA(1,1) processes
% M: nr of observations = 10 000 = nr of rows 
tic;
%              M      N        Method 
garma1=GARMA( 10000, 10, tvec, 1); % Method 1 
toc;
% Elapsed time is 0.016053 seconds. => faster!

tic;
garma2=GARMA( 10000, 10, tvec, 2);% Method 2  
toc;
% Elapsed time is 0.460987 seconds.  => slower 


