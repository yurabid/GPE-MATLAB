% global tcoef omegaBar
% NN0 = task.Ntotal;
% load('TMNdep3.mat');
% tm_coeffs;
% IC = Ja;
% Jb=0;
% invC = Ua*NN0;
R = 1e50;% 380 ~ gamma = 0.005

%%
zavg = zeros(4,200);
niter_outer = 18300;
tfinal = 1830;
times = (0:(niter_outer-1))/(niter_outer-1)*tfinal;
tstep = times(2)-times(1);
IC = 0.0012; 
invC = 1/0.67;%5000*Ua(1); 
averr = zeros(100,100);

%%
 figure;
 hold all;
plot(real(zavg_gpe(5,:))*omegascale/2/pi,real(zavg_gpe(6,:))*omegascale/2/pi);
for i = 12:12
    for ii=[1e-5 22]
        RS = 500*i;
        LR = 10*ii;
        for j=1:200
%     NN = 20*j;
% omegaBar = 1.0;
% omegaPert = 0.05;
% tfinal = (0.06/tcoef + 0.065/2/tcoef)/omegaBar;

% IC = real(interp1(250*(1:30),Ja,NN,'spline','extrap')); % + NN*interp1(250*(1:30),FF,NN,'spline','extrap'));
% invC = NN*interp1(250*(1:30),Ua,NN,'spline','extrap');
% coeff = interp1(250*(1:30),coeffs,NN,'spline','extrap');
% times = (0:niter_outer-1)*ddt*niter_inner;
% tstep = times(2)-times(1);
% coeff = coeffs(1);

% tavg = round(2*pi/(tstep*omegaPert));
% if(tavg>3000)
%     tavg=1;
% end
% zz0=coeff*arrayfun(@rotangle,times,omegaBar+times*0);
% zz0 = 8.6*(2.5e-2/100*j * (1-times/36) .*(times<36) + 0.00*sin(0.2*times));
zz0 = (gather(zavg_gpe(1,j)) * (1-times/36) .*(times<36) + 0.00*sin(0.2*times));
% zz0=coeff*(arrayfun(@rotangle,times,omegaBar+times*0,0.25+times*0) - 0.125*pi);
% zz0=(0.5/invC)*1.01+coeff*0.02*sin(0.5*times);
% zz0=0.005*IC*j*times+0.1*sin(0.02*times);
% zz0=1.36*(1e-9*j)/(rscale/tcoef)*times; % smerzi result

[zztm, phitm] = solve_two_mode(IC,0,invC,R,zz0(1),0,times,zz0,gather(cl(j))/0.98,gather(cr(j))/0.98,RS,LR);

% b = (1/50)*ones(1,50);
% dphitm = filter(b,1,(phitm(2:end)-phitm(1:end-1))/tstep);
% zzz = zz0-zztm;
% plot(zztm);drawnow;
% zavg(j) = sum(zzz(end-tavg+1:end))/tavg;
zavg(:,j) = [zztm(1) zztm(360) zztm(end) sum(zztm(end-5000:end))/5001];
        end
        zint=interp1(zavg(2,:)/0.67,zavg(4,:)/0.67,gather(zavg_gpe(5,:)));
%         averr(i,ii) = sqrt(sum((gather(real(zavg_gpe(6,20:180)))-zint(20:180)).^2));
        plot(zavg(2,:)/0.67*omegascale/2/pi,zavg(4,:)/0.67*omegascale/2/pi); drawnow;
%         averr(i,ii)
    end
end
% save('averr2','averr');