function A = findiff1_arb(r,n,is_periodic)
% build a n-th derivative matrix for arbitrary FD grid r with 5-point stencil
% IN:
%   - r: uniformly growing 1D array of coordinate points
%   - n: order of derivative
%   - is_periodic: type of boundary conditions (
%                  0 - zero Dirichlet (midpoint boundary), 
%                  1 - periodic, 
%                  2 - zero Dirichlet (next point boundary),
%                  3 - zero Neumann (midpoint boundary),
%                  4 - zero Neumann (next point boundary),
%                  5 - continuity condition for singular problems,
% OUT:
%   - A: 2D matrix of the derivative operator

if(nargin==1)
    n=2;
end
if(nargin<=2)
    is_periodic = 0;
end
Ntot = length(r);
r=[4*r(1)-3*r(2) 3*r(1)-2*r(2) 2*r(1)-r(2) r 2*r(end)-r(end-1) 3*r(end)-2*r(end-1) 4*r(end)-3*r(end-1)];
A=zeros(Ntot,Ntot,'like',r);
N_stencil = 5;
max_stencil = (N_stencil-1)/2;
ind = -max_stencil:max_stencil;
for i=4:Ntot+3
    dr = r(i+ind(:))-r(i);
    trans = zeros(N_stencil,N_stencil,'like',r);
    for ii=1:N_stencil
        for jj=1:N_stencil
            trans(ii,jj) = dr(ii)^(jj-1)/factorial(jj-1);
        end
    end
    trans = inv(trans);
    trans = trans(n+1,:);
    tmp = zeros(1,Ntot+6,'like',r); 
    tmp(i+ind(:)) = trans;
    if(is_periodic ==0)
%         tmp(i+ind(:)) = trans;
		tmp(4) = tmp(4) - tmp(3);
        tmp(5) = tmp(5) - tmp(2);
		tmp(Ntot+3) = tmp(Ntot+3) - tmp(Ntot+4);
        tmp(Ntot+2) = tmp(Ntot+2) - tmp(Ntot+5);
%         A(i-2,:) = tmp(3:Ntot+2);
    elseif(is_periodic ==1)
        tmp = zeros(1,Ntot,'like',r);
        tmp(max_stencil+1+ind(:)) = trans;
        tmp = [0,0,0,circshift(tmp,[0,i-max_stencil-4]),0,0,0];
    elseif(is_periodic ==2)
%         tmp = zeros(1,Ntot+4,'like',r);
%         tmp(i+ind(:)) = trans;
		tmp(4) = tmp(4) - tmp(2);
		tmp(Ntot+3) = tmp(Ntot+3) - tmp(Ntot+5);
%         A(i-2,:) = tmp(3:Ntot+2);      
    elseif(is_periodic ==3)
%         tmp = zeros(1,Ntot+4,'like',r);
%         tmp(i+ind(:)) = trans;
		tmp(4) = tmp(4) + tmp(3);
        tmp(5) = tmp(5) + tmp(2);
		tmp(Ntot+3) = tmp(Ntot+3) + tmp(Ntot+4);
        tmp(Ntot+2) = tmp(Ntot+2) + tmp(Ntot+5);
%         A(i-2,:) = tmp(3:Ntot+2);
    elseif(is_periodic ==4)
%         tmp = zeros(1,Ntot+4,'like',r);
%         tmp(i+ind(:)) = trans;
		tmp(4) = tmp(4) + tmp(3) + tmp(2);
		tmp(Ntot+3) = tmp(Ntot+3) + tmp(Ntot+4) + tmp(Ntot+5);
%         A(i-2,:) = tmp(3:Ntot+2);          
	else
%         tmp = zeros(1,Ntot+4,'like',r);
%         tmp(i+ind(:)) = trans;
		tmp(4) = tmp(4) + 2*tmp(3) + 3*tmp(2);
        tmp(5) = tmp(5) - tmp(3) - 2*tmp(2);
		tmp(Ntot+3) = tmp(Ntot+3) + 2*tmp(Ntot+4) + 3*tmp(Ntot+5);
        tmp(Ntot+2) = tmp(Ntot+2) - tmp(Ntot+4) - 2*tmp(Ntot+5);
%         A(i-2,:) = tmp(3:Ntot+2);		   
    end
    A(i-3,:) = tmp(4:Ntot+3);
end
A=sparse(A);
