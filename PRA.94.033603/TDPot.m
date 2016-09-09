% Time-dependent part of the potential
function VVV = TDPot(XXg,YYg,time)

global omegaBar ub1 Ub initphase;

if(nargin == 2)
    addphase = initphase;
    ampl = Ub;
else
    ampl = Ub;
    addphase = rotangle(time,omegaBar) + initphase;
end
cs = cos(addphase+0.5*pi);
sn = sin(addphase+0.5*pi);
% cs2 = cos(-addphase-0.5*pi);
% sn2 = sin(-addphase-0.5*pi);
 tmp = (XXg(end/2+1:end,:)*sn - YYg(end/2+1:end,:)*cs).^2;
VVV = ampl*(exp(-0.1414*tmp));%.*(YYg*sn+XXg*cs>0)+exp(-0.1414*(XXg*sn2 - YYg*cs2).^2).*(YYg*sn2+XXg*cs2>0));
VVV = [flip(VVV,1); VVV];
if(time<36)
    VVV = VVV +  ub1*(1-time/36)*XXg;
%     VVV = VVV +  0.02*sin(0.2*time)*(XXg>0) - 0.02*sin(0.2*time)*(XXg<0);
else
%     VVV = VVV +  0.005*sin(0.2*(time-20))*XXg;
end
end