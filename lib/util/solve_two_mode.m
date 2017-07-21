function [zztm, phitm] = solve_two_mode(IC,Jb,lam,R,z0,phi0,times,edrive)

niter_outer = numel(times);

phitm = zeros(1,niter_outer);
zztm = zeros(1,niter_outer);

zztm(1) = z0;
phitm(1) = phi0;
if(nargin <= 7)
    edrive = times*0;
end
if(numel(edrive) == 1)
    edrive = edrive + times*0;
end
if(R==0)
    G=0;
else
    G=1/R;
end


for i=2:niter_outer
	tstep = times(i)-times(i-1);
%     zztmp = zztm(i-1) - tstep/2*((IC+Jb*zztm(i-1))*sqrt(1-zztm(i-1)^2)*sin(phitm(i-1)) + invC*G*(zztm(i-1)-zz0(i-1)) );
%     phitmp = phitm(i-1) + tstep/2*(invC*(zztm(i-1)-zz0(i-1)) + IC*zztm(i-1)/sqrt(1-zztm(i-1)^2)*cos(phitm(i-1)) + Jb*(2*zztm(i-1)^2-1)/sqrt(1-zztm(i-1)^2)*cos(phitm(i-1)));
%     zztm(i) = zztm(i-1) - tstep*((IC+Jb*zztmp)*sqrt(1-zztmp^2)*sin(phitmp) + invC*G*(zztmp-(zz0(i)+zz0(i-1))/2) ) ;
%     phitm(i) = phitm(i-1) + tstep*(invC*(zztmp-(zz0(i)+zz0(i-1))/2) + IC*zztmp/sqrt(1-zztmp^2)*cos(phitmp) + Jb*(2*zztmp^2-1)/sqrt(1-zztmp^2)*cos(phitmp));
    dmu = (edrive(i-1) + lam*zztm(i-1) + IC*zztm(i-1)/abs(sqrt(1-zztm(i-1)^2))*cos(phitm(i-1)) );
    zztmp = zztm(i-1) - tstep/2*( IC*abs(sqrt(1-zztm(i-1)^2))*sin(phitm(i-1)) + G*dmu);
    phitmp = phitm(i-1) + tstep/2*dmu;
    dmu = ((edrive(i)+edrive(i-1))/2+lam*zztmp + IC*zztmp/abs(sqrt(1-zztmp^2))*cos(phitmp) );
    zztm(i) = zztm(i-1) - tstep*( IC*abs(sqrt(1-zztmp^2))*sin(phitmp)  + G*dmu) ;
    phitm(i) = phitm(i-1) + tstep*dmu;

end

end