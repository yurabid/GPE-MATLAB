if(exist('L','var') ~= 1)
    disp('reinitializing config parameters');
    config
end
r = linspace(-L/2,L/2,N);
rz = linspace(-Lz/2,Lz/2,Nz);
h = r(2)-r(1);
hz = rz(2)-rz(1);
dt = ddt*tcoef*niter_inner;
f=figure;
[XX,YY] = meshgrid(r,r);
angles3 = atan2(XX,YY);

load(sprintf('snapshots/slice_%05d.mat',0));
slice0 = slice(:,:);
%[s0,v0] = velocity_field(slice,r,-(j-1)*dt*0.5*2*pi);
maxx = max(max(abs(slice(:,:))));
[X,Y] = meshgrid(r,r);
%maxx = 39;
for j=1:10000
    if(exist(sprintf('snapshots/slice_%05d.mat',j),'file') ~= 2)
        break
    end
    load(sprintf('snapshots/slice_%05d.mat',j));
    clf(f,'reset')
    time = ddt*niter_inner*j;
%     curangle = (rotangle(time,omegaBar,0.25) + initphase)/pi;
     ang1 = angle(sum(sum(slice.*(XX>0))));% angle(sum(sum(slice.*(angles3>(-pi + curangle*pi)).*(angles3<-curangle*pi))));
     slice = slice.*exp(-1i*ang1);
	[coresp, coresm] = detect_core(slice,X,Y);
% 	[s,v] = velocity_field(slice,r,-(j-1)*dt*1.5*2*pi,slice0);
%     [s,v] = velocity_field(slice,r,-angles(j)*pi-pi/2,slice0);
	axes('Position',[0.2,0.1,0.62,0.84])
    imagesc(r,r,angle(slice));
    set(gca,'YDir','normal');
	set(f, 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 6]);
    annotation(f,'textbox',...
        [0.1246428571428571 0.460428571428572 0.0482142857142857 0.0833333333333334],...
        'String',{'Y'},...
        'FontSize',14,...
        'FitBoxToText','on',...
        'LineStyle','none');
    annotation(f,'textbox',...
        [0.496785714285713 0 0.0482142857142857 0.0690476190476191],...
        'String',{'X'},...
        'FontSize',14,...
        'FitBoxToText','on',...
        'LineStyle','none');
    annotation(f,'textbox',...
        [0.6 0.788961038961044 0.21092023633678 0.127424242424244],...
        'String',{sprintf('t = %0.4f',(j*dt))},...
        'FontWeight','bold',...
        'FontSize',20,...
        'FitBoxToText','on',...
        'LineStyle','none',...
        'Color',[0 0 0]);
	hold on;
% 	s2 = streamline(s);
% 	set(s2,'Color',[1,1,1]);
% 	set(s2,'LineWidth',1);
	plot(coresp(:,1),coresp(:,2),'MarkerSize',6,'Marker','+','LineStyle','none','Color',[1 0 0]);
    plot(coresm(:,1),coresm(:,2),'MarkerSize',6,'Marker','o','LineStyle','none','Color',[0 0 1]);
    print(f,'-dpng',sprintf('snapshots/phase_%05d.png',j), '-r150');
    image(r,r,abs(slice)/maxx*64);
    set(gca,'YDir','normal');
    annotation(f,'textbox',...
        [0.6 0.788961038961044 0.21092023633678 0.127424242424244],...
        'String',{sprintf('t = %0.4f',(j*dt))},...
        'FontWeight','bold',...
        'FontSize',20,...
        'FitBoxToText','on',...
        'LineStyle','none',...
        'Color',[1 1 0]);
	hold on;
% 	s2 = streamline(s);
% 	set(s2,'Color',[1,1,1]);
% 	set(s2,'LineWidth',1);
	plot(coresp(:,1),coresp(:,2),'MarkerSize',6,'Marker','+','LineStyle','none','Color',[1 0 0]);
    plot(coresm(:,1),coresm(:,2),'MarkerSize',6,'Marker','o','LineStyle','none','Color',[0 0 1]);
	print(f,'-dpng',sprintf('snapshots/sol_%05d.png',j), '-r150');
end