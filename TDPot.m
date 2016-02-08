% Time-dependent part of the potential
function VVV = TDPot(time)

global XXg Ub vv linit rscale;

% vv = 0.004;
ampl = Ub;
VVV = ampl*exp(-2*((XXg - vv*time - linit)./(3.5e-6/rscale)).^2);

end