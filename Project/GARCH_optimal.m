function [gbest,abest,ICmat ] = GARCH_optimal( data, gmax, amax, arp, map, method)
%This function computes the optimal AR and MA lags for an ARMA(p,q) model
%based on information criteria. Inputs:
% data: T-by-1 vector of data to be modeled by an ARMA(p,q) model (log ret)
% gmax: integer. Maximum GARCH lags considered in the model
% amax: integer. Maximum ARCH lags considered in the model
% arpar: number of AR parameters in the mean equation
% mapar: number of MA parameters in the mean equation
% method: string. Methods to used for choosing the optimal model. Choice:
%'AIC': Akaike Information Criterion
%'BIC': Bayesian Information Criterion
%'HQIC': Hannan-Quinn Information Criterion
%Outputs:
%gbest: optimal GARCH lags
%abest: optimal ARCH lags.
%ICmat: Matrix of the chosen information criteria



T=length(data);
%Restrict the possible number of parameters:
% amax = maximum ARCH lags considered in the model
% gmax = maximum GARCH lags considered in the model
if amax+1+gmax>=T/10 % = 2517/10 = 251.7
    % Here amax + 1 g gmax = 2+1+2 = 5 < 251.7
    error('Too many parameters. Reduce pmax and qmax');
end

ICmat=zeros(1+gmax,amax); % (2 + 1) x 2 = 3 x 2 matrix

%start from GARCH(0,0), which is an ARMA(arp, arq) model). Note that at
%least an ARCH term is required for a GARCH model.
for i=1:gmax+1 % i = 1,2,3
    for j=1:amax+1 % j = 1,2,3
        %If it is GARCH(0,0), treat it as an ARMA model
        if i*j==1 % if i = 1, and j = 1
                  % equivalent with GARCH(0,0)
            Mdltemp=arima('ARLags',arp, 'MAlags', map); % ARIMA(arp,map)
            [EstMdl,EstParamCov,logL,info] = estimate(Mdltemp,data,'Display','off');
        elseif i>1 
            if j>1 % if not a GARCH(0,0) model = ARIMA(arp,map)
                Mdltemp=arima('ARLags',arp, 'MAlags', map,'Variance',garch(i-1,j-1));
                [EstMdl,EstParamCov,logL,info] = estimate(Mdltemp,data,'Display','off');
            else
                %If it is a GARCH(i,0) model, return negative infinity as
                %likelihood
                logL=-Inf;
            end
        end
                
        if strcmp(method,'AIC')
            ICmat(i,j)=-2*logL/T+2*(arp+map+i+j)/T;
        elseif strcmp(method,'BIC')
            ICmat(i,j)=-2*logL/T+(arp+map+i+j)*log(T)/T;
        elseif strcmp(method,'HQIC')
            ICmat(i,j)=-2*logL/T+2*(arp+map+i+j)*log(log(T))/T;
        else
            error('Unknown choice of information criterion');
        end
    end
end
%Find the location of the minimum information
[gbest,abest]=find(ICmat==min(min(ICmat)));
gbest=gbest-1;
abest=abest-1;

end

