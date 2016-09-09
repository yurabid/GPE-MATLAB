% figure;
% hold all;
muavg = zeros(1,100);
    omp = 0.04;
    IC = 0.0012;
    C = 0.67;
    LP = 0;
    RS = 1;
    z0 = @(t,u0) u0*(1-t/20).*(t<20);
for i = 1:100
    u0 = i/100;
    [t,y] = ode45(@(t,y) [y(4); -IC*sqrt(1-y(2).^2).*sin(y(3)); -LP*y(4) + (y(2)-z0(t,u0))/C; -y(1)*omp - RS*omp*y(4) + IC*omp*sin(y(3))], (0:0.1:2048), [0; u0; 0; 0]);
%     [t,y] = ode45(@(t,y) [y(2); IC/C*sin(y(1))-y(3)/C; y(2)/LP - RS/LP*y(3) - omp*y(4); y(3)], (0:3000), [0; u0; 0; 0]);
    mmu = y(:,2)/C-y(:,4)*LP;
    muavg(i) = sum(mmu(end-2000:end))/2001;
%     plot(y(:,1));
end
% plot(muavg)