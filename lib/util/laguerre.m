% y=Laguerre(n,m,x) returns the (generalised) Laguerre Polynomial
% L_{n}^{m}(x) by calculating the coefficients of the Taylor series
% following http://en.wikipedia.org/wiki/Laguerre_polynomials and then
% using polyval to calculate the actual function values. Required are input
% n>0.
% The idea is based upon the solution in
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=4914
% however without distinguishing m=0 from other polynomials.

function y=laguerre(n,m,x)
P = zeros(n+1,1);
for v=0:n
P(n+1-v) = (-1)^v * factorial(n+m)/factorial(n-v)/factorial(m+v)/factorial(v);
end;
y = polyval(P,x);
