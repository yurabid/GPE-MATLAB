function A = findiff1_arb(r,n,is_periodic)
% build a second derivative matrix for arbitrary FD grid r with 3-point stencil
if(nargin==1)
    n=2;
end
if(nargin<=2)
    is_periodic = 0;
end
Ntot = length(r);
r=[3*r(1)-2*r(2) 2*r(1)-r(2) r 2*r(end)-r(end-1) 3*r(end)-2*r(end-1)];
A=zeros(Ntot,Ntot);
N_stencil = 3;
max_stencil = (N_stencil-1)/2;
ind = -max_stencil:max_stencil;
for i=3:Ntot+2
    dr = r(i+ind(:))-r(i);
    trans = zeros(N_stencil,N_stencil);
    for ii=1:N_stencil
        for jj=1:N_stencil
            trans(ii,jj) = dr(ii)^(jj-1)/factorial(jj-1);
        end
    end
    trans = inv(trans);
    trans = trans(n+1,:);
    if(is_periodic ==0)
        tmp = zeros(1,Ntot+4);
        tmp(i+ind(:)) = trans;
        A(i-2,:) = tmp(3:Ntot+2);
    else
        tmp = zeros(1,Ntot);
        tmp(max_stencil+1+ind(:)) = trans;
        A(i-2,:) = circshift(tmp,[0,i-max_stencil-3]);   
    end
end
A=sparse(A);
