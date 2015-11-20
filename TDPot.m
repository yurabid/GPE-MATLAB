% Time-dependent part of the potential
function VVV = TDPot(time)

global tcoef XXg YYg omegaBar ub1 Ub addphase;

% Parameters of time-dependent potential
% ub1 = 0.53 * 2.6235; % MU = 2.6235 for 2.5e4 particles
% ub1 = 0.53 * 11.7552; % MU = 11.7552 for 7.5e5 particles
% tz = 1 / tcoef;
omegaBar = 1; %7.26467;
t0 = 0; %0.1/tcoef;
t1 = t0+0.06/(omegaBar)/tcoef;
t2 = t1+0.005/(omegaBar*2)/tcoef;
t3 = t2+0.06/(omegaBar)/tcoef;
beta = omegaBar*2*pi*tcoef/(t1-t0);
ampl = 0;

if(time < t0)
    ampl = Ub/2 + Ub/2*time/t0;
    addphase = 0;
elseif(time < t1)
    ampl = Ub;
    addphase = beta*(time-t0)^2/2;
elseif(time<t2)
    ampl = Ub;
    addphase = 0.06*pi + omegaBar*2*pi*tcoef*(time-t1);
elseif(time<t3)
    ampl = Ub;
    addphase = 0.065*pi + omegaBar*2*pi*tcoef*(time-t2) - beta*(time-t2)^2/2; 
    
else
    ampl = Ub;
    addphase = 0.125*pi;
%     ww = omegaBar*2*pi*tcoef;
%     addphase = pi/90*sin(3.5*2*pi*tcoef*time)+ww*(time); %-0.2*tz);
end
% if(addphase > 0.5*pi)
%     addphase = 0.5*pi;
% end

% if(time < 0.2*tz)
%     coef = ub1*time/(0.2*tz);
% else
%     coef = ub1;
% end
% addphase = 0.125*pi;
cs = cos(addphase+0.5*pi);
sn = sin(addphase+0.5*pi);
cs2 = cos(-addphase-0.5*pi);
sn2 = sin(-addphase-0.5*pi);
VVV = ampl*(exp(-0.1414*(XXg*sn - YYg*cs).^2).*(YYg*sn+XXg*cs>0)+exp(-0.1414*(XXg*sn2 - YYg*cs2).^2).*(YYg*sn2+XXg*cs2>0));
end