function [pbest,qbest ] = ARMA_optimal( data, pmax, qmax, method )
%This function computes the optimal AR and MA lags for an ARMA(p,q) model
%based on information criteria. Inputs:
% data: T-by-1 vector of data to be modeled by an ARMA(p,q) model
% pmax: integer. Maximum AR lags considered in the model
% qmax: integer. Maximum MA lags considered in the model
% method: string. Methods to used for choosing the optimal model. Choice:
%'AIC': Akaike Information Criterion
%'BIC': Bayesian Information Criterion
%'HQIC': Hannan-Quinn Information Criterion
%Outputs:
%pbest: optimal AR lags
%qbest: optimal MA lags.

T=length(data);
%Restrict the possible number of parameters:
if pmax+1+qmax>=T/10
    error('Too many parameters. Reduce pmax and qmax'); 
end

ICmat=zeros(pmax+1,qmax+1);


for i=1:pmax+1
    for j=1:qmax+1
        Mdltemp=arima(i-1,0,j-1);
       [~,~,logL,~] = estimate(Mdltemp,data,'Display','off'); 
       if strcmp(method,'AIC')
       ICmat(i,j)=-2*logL/T+2*(-1+i+j)/T;
       elseif strcmp(method,'BIC')
       ICmat(i,j)=-2*logL/T+(-1+i+j)*log(T)/T;
       elseif strcmp(method,'HQIC')
       ICmat(i,j)=-2*logL/T+2*(-1+i+j)*log(log(T))/T;
       else
           error('Unknown choice of information criterion');
       end
    end
end
%Find the location of the minimum information
[pbest,qbest]=find(ICmat==min(min(ICmat)));
pbest=pbest-1;
qbest=qbest-1;

end

