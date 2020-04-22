function [ret]=myskewness(x)
ret=mean((x-mean(x)).^3)./var(x).^(3/2);
end