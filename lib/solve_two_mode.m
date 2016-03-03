function [zztm, phitm] = solve_two_mode(IC,Jb,invC,R,z0,phi0,times,zz0)

niter_outer = numel(times);

phitm = zeros(1,niter_outer);
zztm = zeros(1,niter_outer);

zztm(1) = z0;
phitm(1) = phi0;

for i=2:niter_outer
	tstep = times(i)-times(i-1);
    zztmp = zztm(i-1) - tstep/2*((IC+Jb*zztm(i-1))*sqrt(1-zztm(i-1)^2)*sin(phitm(i-1)) + invC/R*(zztm(i-1)-zz0(i-1)) );
    phitmp = phitm(i-1) + tstep/2*(invC*(zztm(i-1)-zz0(i-1)) + IC*zztm(i-1)/sqrt(1-zztm(i-1)^2)*cos(phitm(i-1)) + Jb*(2*zztm(i-1)^2-1)/sqrt(1-zztm(i-1)^2)*cos(phitm(i-1)));
    zztm(i) = zztm(i-1) - tstep*((IC+Jb*zztmp)*sqrt(1-zztmp^2)*sin(phitmp) + invC/R*(zztmp-(zz0(i)+zz0(i-1))/2) ) ;
    phitm(i) = phitm(i-1) + tstep*(invC*(zztmp-(zz0(i)+zz0(i-1))/2) + IC*zztmp/sqrt(1-zztmp^2)*cos(phitmp) + Jb*(2*zztmp^2-1)/sqrt(1-zztmp^2)*cos(phitmp));
end

end