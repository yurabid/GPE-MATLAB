function [y, errors] = lerch(z,s,a) 
%% polylog - Computes the Lerch transcendent: \Phi(z,s,a)

if nargin~=3
    errors=1;
    error('[Error in: lerch function] Inappropriate number of input arguments!')
end

% display more digits in Matlab terminal:
%format long

nsum = 1000;
y = 0.*z;
tmp = 1;

for i=0:nsum
    y = y + tmp./(i+a).^s;
    tmp=tmp.*z;
end

