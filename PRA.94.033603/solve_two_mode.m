function [zztm, phitm] = solve_two_mode(IC,Jb,invC,R,z0,phi0,times,zz0,cl,cr,RS,LP1)

niter_outer = numel(times);

phitm = zeros(1,niter_outer);
zztm = zeros(1,niter_outer);
up = zeros(1,niter_outer);
ip = zeros(1,niter_outer);
up2 = zeros(1,niter_outer);
ip2 = zeros(1,niter_outer);

%     RS = 3000;
    omp = cl^2;%0.23^2;
    omp2 = cr^2;
%     LP1 = 150;
    LP2 = LP1;%*omp/omp2;
    CP1 = 1.0/(LP1*omp);
    CP2 = 1.0/(LP2*omp2);
zztm(1) = z0;
phitm(1) = phi0;

for i=2:niter_outer
	tstep = times(i)-times(i-1);
    zztmp = zztm(i-1) - tstep/2*((IC+Jb*zztm(i-1))*sqrt(1-zztm(i-1)^2)*sin(phitm(i-1)) + invC/R*(zztm(i-1)-zz0(i-1)) );
    phitmp = phitm(i-1) + tstep/2*(-ip(i-1)-ip2(i-1)+invC*(zztm(i-1)-zz0(i-1)) + IC*zztm(i-1)/sqrt(1-zztm(i-1)^2)*cos(phitm(i-1)) + Jb*(2*zztm(i-1)^2-1)/sqrt(1-zztm(i-1)^2)*cos(phitm(i-1)));
    uptmp = up(i-1) + tstep/2*(ip(i-1));
    iptmp = ip(i-1) + tstep/2*( (IC+Jb*zztm(i-1))*sqrt(1-zztm(i-1)^2)*sin(phitm(i-1)) - ip(i-1)/(RS) - up(i-1)/LP1)/CP1;
    up2tmp = up2(i-1) + tstep/2*(ip2(i-1));
    ip2tmp = ip2(i-1) + tstep/2*( (IC+Jb*zztm(i-1))*sqrt(1-zztm(i-1)^2)*sin(phitm(i-1)) - ip2(i-1)/(RS) - up2(i-1)/LP2)/CP2;
    
    zztm(i) = zztm(i-1) - tstep*((IC+Jb*zztmp)*sqrt(1-zztmp^2)*sin(phitmp) + invC/R*(zztmp-(zz0(i)+zz0(i-1))/2) ) ;
    phitm(i) = phitm(i-1) + tstep*(-iptmp-ip2tmp+invC*(zztmp-(zz0(i)+zz0(i-1))/2) + IC*zztmp/sqrt(1-zztmp^2)*cos(phitmp) + Jb*(2*zztmp^2-1)/sqrt(1-zztmp^2)*cos(phitmp));
    up(i) = up(i-1) + tstep*(iptmp);
    ip(i) = ip(i-1) + tstep*( (IC+Jb*zztmp)*sqrt(1-zztmp^2)*sin(phitmp) - iptmp/(RS) - uptmp/LP1)/CP1;
    up2(i) = up2(i-1) + tstep*(ip2tmp);
    ip2(i) = ip2(i-1) + tstep*( (IC+Jb*zztmp)*sqrt(1-zztmp^2)*sin(phitmp) - ip2tmp/(RS) - up2tmp/LP2)/CP2;
end
% zztm = zztm - LP*ip/invC - LP*ip2/invC;
% plot([z0/0.67/cl],[sum(abs(ip)+abs(ip2))],'+');drawnow;
end