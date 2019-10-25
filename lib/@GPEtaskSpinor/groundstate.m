function [phi, varargout] = groundstate(task,eps,phi)
% groundstate - Calculate the stationary state of multicomponent GPE with split-step Imaginary Time Propagation method.
%
%  Usage :
%    phi = task.groundstate_itp(eps)
%    phi = task.groundstate_itp(eps,phi0)
%    [phi, mu] = task.groundstate_itp(eps)
%    [phi, mu, mu2] = task.groundstate_itp(eps)
%    [phi, mu, mu2, E] = task.groundstate_itp(eps)
%  Input
%    eps   :  desired accuracy (applied to chemical potential)
%    phi0  :  initial approximation of the wave function
%  Output
%    phi      :  calculated stationary state
%    mu       :  array of chemical potential values from renormalization
%    mu2      :  array of chemical potential values from integral evaluation
%    E        :  array of energy values

grid = task.grid;
% rscale = (grid.x(2)-grid.x(1));
rscale = grid.weight^(1/grid.ndims);
escale = 1/rscale^2;
tscale = rscale^2;
phiscale = rscale^(-grid.ndims/2);
kk = grid.kk/escale;
ncomp = size(task.g,1);
task.ncomp = ncomp;
VV = task.getVtotal(0);
V = VV{1};
dt = 2;
coupl = task.coupling/escale;
g = task.g;
if(nargin <= 2)
    phi = cellfun(@(x) grid.normalize(rand(size(V),'like',V) + 1i*rand(size(V),'like',V))./sqrt(ncomp/task.Ntotal)/phiscale, cell(1,ncomp), 'UniformOutput', false);
else
    ntot = real(grid.integrate(phi{1}.*conj(phi{1})) + grid.integrate(phi{2}.*conj(phi{2})));
    phi{1} = phi{1}/phiscale*sqrt(task.Ntotal/ntot);
    phi{2} = phi{2}/phiscale*sqrt(task.Ntotal/ntot);
end
ekk = exp(-kk*dt*0.5);
MU = zeros(1000,1,'like',V);
MU2 = zeros(1000,1,'like',V);
dts = zeros(1000,1,'like',V);
EE = zeros(1000,1,'like',V);
i = 1;
iswitch = 50;

tmp2 = cell(1,ncomp);
cang = angle(coupl);
cosom = cosh(dt*abs(coupl));
sinomm = sinh(dt*abs(coupl)).*exp(-1i*cang);
sinomp = sinh(dt*abs(coupl)).*exp(1i*cang);

while true
    for j =1:ncomp
        tmp2{j} = VV{j}/escale;
        for k = 1:ncomp
            tmp2{j} = tmp2{j} + g(j,k)*abs(phi{k}).^2;
        end
    end
    if(task.nlincpl~=0)
        tmp3 = coupl + task.nlincpl*(abs(phi{1}).^2+abs(phi{2}).^2);
        cosom = cosh(dt*tmp3);
        sinomm = sinh(dt*tmp3);
        sinomp = sinomm;
    end
    phi{1} = grid.ifft(ekk.*grid.fft(phi{1}));
    phi{2} = grid.ifft(ekk.*grid.fft(phi{2}));
    tmp2{1} = exp(-tmp2{1}*dt).*(cosom.*phi{1} - sinomp.*phi{2});
    tmp2{2} = exp(-tmp2{2}*dt).*(cosom.*phi{2} - sinomm.*phi{1});
    phi{1} = grid.ifft(ekk.*grid.fft(tmp2{1}));
    phi{2} = grid.ifft(ekk.*grid.fft(tmp2{2}));
    
    ntot = real(grid.integrate(phi{1}.*conj(phi{1}) + phi{2}.*conj(phi{2})))*phiscale^2;
    mu = sqrt(task.Ntotal/ntot);
    dts(i) = dt*tscale;
    MU(i) = log(mu)/dts(i);
    phi=task.mul(phi,mu);
    
    if(nargout >= 3)
        MU2(i) = real(task.inner(task.mul(phi,phiscale),task.applyham(task.mul(phi,phiscale))))/task.Ntotal;
    else
        MU2(i) = MU(i);
    end
    if(nargout >= 4)
        EE(i) = task.get_energy(task.mul(phi,phiscale))/task.Ntotal;
    else
        EE(i) = MU2(i);
    end
    if((i-iswitch)>200 && mod(i,10) == 0)
        delta = abs((EE(i)-EE(i-10))^2/(EE(i)-2*EE(i-10)+EE(i-20)))/EE(i);
%         delta = (abs(EE(i)-EE(i-20)))/dt/20;
        if(delta<dt*0.0000001 || delta < eps)
            if (dt<eps || dt<2e-3)
                break;
            else
                dt = dt/1.5;
                ekk = exp(-kk*dt*0.5);
                cosom = cosh(dt*abs(coupl));
                sinomm = sinh(dt*abs(coupl)).*exp(-1i*cang);
                sinomp = sinh(dt*abs(coupl)).*exp(1i*cang);
                iswitch = i;
            end
        end
    end
    if(i>=50000)
        warning('Convergence not reached');
        break;
    end
    i=i+1;
end

if(nargout >= 2)
    MU = MU(1:i-1);
    varargout{1} = MU;
end
if(nargout >= 3)
    MU2 = MU2(1:i-1);
    varargout{2} = MU2;
end
if(nargout >= 4)
    EE = EE(1:i-1);
    varargout{3} = EE;
end
phi = task.mul(phi,phiscale);
task.init_state = phi;
task.current_state = phi;
end
