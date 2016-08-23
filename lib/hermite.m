function H = hermite(n,x,nflag)
%HERMITE Hermite polynomials of order N.
%   H = HERMITE(N,X) returns the Hermite polynomials of orders N at
%   locations X. N must be a vector of positive integers, and X must be a
%   vector of real numbers.
%
%   H = HERMITE(N,X,'norm') returns the normalized form of the Hermite 
%   polynomials, with normalization such that
%
%     integral{ x=(-inf,inf); exp(-x^2) [ H_n(x) ]^2 } dx = 1
% 
%   Example:
% 
%       % Plot the first five Hermite polynomials
%       x = -3:0.05:3;
%       n = 0:4;
%       y = hermite(n,x);
%       plot(x,y,'LineWidth',2)
%       xlabel('X')
%       ylabel('H_n(X)')
%       title('Hermite polynomials for N = 0 to 4')
%       set(gca,'YLim',[-25 25])
%       legend(strcat('N=',num2str(n')),'Location','SouthEast')

% Paul Fricker 5/25/2012
% Copyright 2012 MathWorks

% Check the inputs:
% -----------------
if ( ~any(size(n)==1) )
    error('hermite:NVector','N must be a vector (or scalar).')
end
if any(n<0)
    error('hermite:NPositive','All elements of N must be positive.')
end
if any(round(n)~=n)
    error('hermite:NIntegers','All elements of N must be integers.')
end

if ~any(size(x)==1)
    error('hermite:XVector','X must be a vector (or scalar).')
end
if any(~isreal(x))
    error('hermite:XReal','All elements of X must be real.')
end

% Check normalization:
% --------------------
if nargin==3
    if ~ischar(nflag)
        error('hermite:NormFlag', ...
              'The normalization flag must be the character string ''norm''.')
    else
        isnorm = strcmpi(nflag,'norm');
        if ~isnorm
            error('hermite:Normalization','Unrecognized normalization flag.')
        end
    end
else
    isnorm = false;
end

% Initialize a few variables:
% ---------------------------
n = n(:);
fn = floor(n/2);
x = 2*x(:);
length_x = length(x);
length_n = length(n);

% Pre-compute the required powers of x:
% -------------------------------------
p = arrayfun(@(N,M)(N - 2*(0:M)),n,fn,'UniformOutput',false);
p = unique([p{:}]);

% Pre-compute the values of x raised to the required powers, and
% compile them in a matrix:
% -------------------------
if p(1)==0 % Don't compute x.^0
    xp = arrayfun(@(P)x.^P,p(2:end),'UniformOutput',false);
    xp = [ones(length_x,1,'like',x) xp{:}];
else
    xp = arrayfun(@(P)x.^P,p,'UniformOutput',false);
    xp = [xp{:}];
end

% Compute the Hermite polynomials:
% --------------------------------
H = zeros(length_x,length_n,'like',x);
for k = 1:length(n)
    for m = 0:fn(k)
        is_the_power = p == (n(k) - 2*m);
        H(:,k) = H(:,k) + (1 - 2*mod(m,2))/ ...      % (-1)^m
                          prod(2:m)/ ...             % factorial(m)
                          prod(2:(n(k) - 2*m))* ...  % factorial(n-2m)
                          xp(:,is_the_power);        % (2x)^(n-2m)
    end
    H(:,k) = prod(2:n(k))*H(:,k);
    
    if isnorm
        H(:,k) = H(:,k)/sqrt(prod(2:n(k)) * 2^n(k) * sqrt(pi));
    end
end

% EOF hermite