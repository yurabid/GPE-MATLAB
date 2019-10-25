% Time-dependent part of the potential
function VVV = TDPot(XXg,~,time)

global Ub vv rscale;

VVV = Ub*exp(-2*((XXg - vv*time)./(3.5e-6/rscale)).^2); % Gaussian barrier with amplitude Ub and 1/e^2 half width 3.5 \mu m

end