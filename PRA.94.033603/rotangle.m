function addphase = rotangle(time,omegaBar, thetatotal)
global tcoef;
thetaacc = 0.06;
thetadecc = 0.00;
omegaPert = 0.00;
if(omegaBar==0)
    addphase=0;
else
    if(nargin == 3)
        thetaconst = thetatotal - 0.06;
        if(thetaconst<0)
            thetaacc = thetatotal/2;
            thetadecc = thetatotal/2;
            thetaconst = 0;
        end
    else
        thetaconst = 0.065;
    end
    t0 = 0; %0.1/tcoef;
    t1 = t0+thetaacc/(omegaBar+1e-100)/tcoef;
    t2 = t1+thetaconst/(omegaBar*2+1e-100)/tcoef;
    t3 = t2+thetadecc/(omegaBar+1e-100)/tcoef;
    beta = omegaBar*2*pi*tcoef/(t1-t0);
    %     ampl = 0;
    
    if(time < t0)
        addphase = 0;
    elseif(time < t1)
        addphase = beta*(time-t0)^2/2;
    elseif(time<t2)
        addphase = thetaacc*pi + omegaBar*2*pi*tcoef*(time-t1); % + 0.05*sin(12*2*pi*tcoef*(time-t1));
    elseif(time<t3)
        %     ampl = Ub;
        addphase = (thetaacc+thetaconst)*pi + omegaBar*2*pi*tcoef*(time-t2) - beta*(time-t2)^2/2;
        %
    else
        %     ampl = Ub;
        addphase = (thetaacc + thetaconst + thetadecc)*pi; % + pi/200*sin(89*2*pi*tcoef*(time-t1));
        
        %     ww = omegaBar*2*pi*tcoef;
        %     addphase = pi/90*sin(3.5*2*pi*tcoef*time)+ww*(time); %-0.2*tz);
    end
end
addphase = addphase  + 0.005*sin(omegaPert*time);
% if(addphase > 0.5*pi)
%     addphase = 0.5*pi;
% end

% if(time < 0.2*tz)
%     coef = ub1*time/(0.2*tz);
% else
%     coef = ub1;
% end
% addphase = 0.125*pi;
end
