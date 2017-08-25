function [Nc,mutot,ecut] = sgpe_params(omegabar,a,Tdim,N)

eps = 1e-10;
omz = omegabar.^3;
e0=1+omz/2;
abar = 1/sqrt(omegabar);
ecuthf = Tdim*log(1+1);
ecgrid = (1:1000)/1000*ecuthf;
ecgrid2 = (1:100000)/10;
dec = ecgrid(2)-ecgrid(1);
dec2 = ecgrid2(2)-ecgrid2(1);
delta = 1;
Nr = 1.5*N;
Nl = 1;
i=0;
while delta>eps
    Nc = (Nr+Nl)/2;
    delta = (Nr-Nl)/N;
    mutot = omegabar/2*(15*Nc*a/abar).^(2/5);
    
%     ecut = Tdim*log(1+1/3) + mutot; %thermal population of 2 atoms per mode at Ecut %26.35*omegabar;
%     Nnc = polylog_inc(3,exp(mutot./Tdim),ecut./Tdim)./(omegabar./Tdim).^3;

    I = sum(dens_states(ecgrid,mutot,omegabar)).*dec;
    ecut = (I*6*omegabar^3 + e0^3)^(1/3);
    Nnc = sum(abs(dens_states(ecgrid2,mutot,omegabar))./(exp(ecgrid2/Tdim)-1)).*dec2;
    
    Nest =  Nnc+Nc;
    if(Nest>N)
        Nr=Nc;
    else
        Nl=Nc;
    end
end