function ret = HAR_eval(arg)
%This function computes the MAPE, MSPE and QLIKE of the forecasts from
%different series. 
%Input: arg is a (T-by-2*C) matrix of the following form:
%[true forecast1 true forecast2 true forecast3 ... true forecastC ]
%Output: ret is a (C-by-3) matrix:
%[MAPE1 MSPE1 QLIKE1; MAPE2 MSPE2 QLIKE2 ; ... ; MAPEC MSPEC QLIKEC]

%Get the number of columns in the arg
C=size(arg,2);
%Set an index for the forecasts that retrieves all the even columns from
%arg
ind=2:2:C;
%Initialize the return series. 
ret=zeros(C/2,3);
%Computing MAPE
ret(:,1)=mean(abs(arg(:,ind)-arg(:,ind-1)))';
%Computing MSPE
ret(:,2)=mean((arg(:,ind)-arg(:,ind-1)).^2)';
%Computing QLIKE
ret(:,3)=mean(arg(:,ind)./arg(:,ind-1)-log(arg(:,ind)./arg(:,ind-1))-1)';
end

