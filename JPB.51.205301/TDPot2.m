% Time-dependent part of the potential
function VVV = TDPot(XXg,YYg,time)

global Ub vv linit rscale ub1 tshift;

% vv = 0.004;
ampl = Ub;
%VVV = ampl*exp(-2*((XXg - vv*time - linit)./(3.5e-6/rscale)).^2); %smerzi
VVV = ampl*exp(-2*((XXg)./(0.7e-6/rscale)).^2); %levi
if(time==0)
%     VVV = VVV + 1*(XXg>0);
else
%     VVV = VVV + 0.02*sin(0.2*time)*(XXg>0);
end
if(time<tshift)
    VVV = VVV +  ub1*(1-time/tshift)*XXg;
% else
%     VVV = VVV +  0.03*sin(0.5*(time-50))*XXg;
end
end