%This function computes the log difference between the i-th row and the
%(i-1)-th row of a matrix for i=2:end
%The first element is padded with zeroes

function ret=dlog(X) % returnerer ret - så ret er det som returneres 
                     % dlog = navnet på funksjonen
                     % X: argumentet
c=size(X,2); % data matrise har 2 dimensjoner; den først er antall rader som er 7742
             % den andre er antall kolonner som er 9 = c
ret=[zeros(1,c);diff(log(X))]; % zeros(1,c)=zeros(1,9) blir en 1x9 vektor som returneres
                               % 
end                  %  
% log returns = log(p_t - p_(t-1)) = log(p_t) - log(p_t-1)
%diff(X) og X er en vektor av lengde m => Y=diff(X)
%Y = [X(2)-X(1) X(3)-X(2) ... X(m)-X(m-1)] - 1x(m-1) matrise 
% X: pxm matrise => 
%  Y = diff(X)= [X(2,:)-X(1,:); X(3,:)-X(2,:); ... X(p,:)-X(p-1,:)] 
% elements are the differences between the rows of X.
