function h = Hermite(n,x)

h = zeros(1,n+1);

n_fact = factorial(n);

for m=0:floor(n/2)
    h(2*m+1) = n_fact * (-1)^m / (factorial(m) * factorial(n-2*m)) * 2^(n-2*m);
end

if exist('x','var')
    h = polyval (h, x);
end