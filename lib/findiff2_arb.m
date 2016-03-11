function A = findiff2_arb(r,n)
% build a second derivative matrix for arbitrary FD grid r with 5-point stencil
if(nargin==1)
    n=2;
end
r=[3*r(1)-2*r(2) 2*r(1)-r(2) r 2*r(end)-r(end-1) 3*r(end)-2*r(end-1)];
Ntot = length(r);
A=zeros(Ntot,Ntot);
ind = [-2 -1 0 1 2];
for i=3:Ntot-2
    dr = r(i+ind(:))-r(i);
    N_stencil = length(ind);
    trans = zeros(N_stencil,N_stencil);
    for ii=1:N_stencil
        for jj=1:N_stencil
            trans(ii,jj) = dr(ii)^(jj-1)/factorial(jj-1);
        end
    end
    trans = inv(trans);
    trans = trans(n+1,:);
    tmp = zeros(1,Ntot);
    tmp(i+ind(:)) = trans;
    A(i,:) = tmp;   
end
A=sparse(A(3:end-2,3:end-2));
