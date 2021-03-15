function [y, errors] = polylog_inc(n,z,x) 
%% polylog - Computes the incomplete polylogarithm of z: Li_n(z)

if nargin~=3
    errors=1;
    error('[Error in: polylog_inc function] Inappropriate number of input arguments!')
end

% display more digits in Matlab terminal:
%format long

nsum = 100;
y = zeros(size(z),'like',z);
tmp = ones(size(z),'like',z);

for i=1:nsum
    tmp=tmp.*z; 
    y = y + tmp./i.^n.*gammainc(x.*i,n,'upper');
end

