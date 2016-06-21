if(exist('L','var') ~= 1)
    disp('reinitializing config parameters');
    config
end
r = linspace(-L/2,L/2,N);
rz = linspace(-Lz/2,Lz/2,Nz);
h = r(2)-r(1);
hz = rz(2)-rz(1);
dt = ddt*tcoef*niter_inner;
f = figure;
%dt = 0.004*10/(100*2*pi);
% VG = gather(V(:,:,16));
% for i = [3 111 124 343 3] - Fig. 5
% for i = [3 430 503 999 3] - Fig. 6
% for i = [64 440 552 900 64] - Fig. 8
% for i = [3 4009 3] - Fig. 7
for i = [1 55 70 471 1]
    load(sprintf('snapshots_m_1_s_2_om_0_25/core_%05d.mat',i));
    clf(f,'reset');
    axes1 = axes('Parent',f,'FontSize',24,...
    'Position',[0.129224021592444 0.0701174947216677 0.850701754385965 0.914101886166624]);
    xlim(axes1,[-9 9]);
    ylim(axes1,[-9 9]);
    zlim(axes1,[-7 7]);
view(axes1,[40.5 24]);
%view(axes1,[0 0]);

    box(axes1,'on');
    grid(axes1,'off');
    hold(axes1,'all');
    %plot3(cores(:,1),cores(:,2),cores(:,3),'MarkerSize',2,'Marker','o','LineStyle','none');
%     maskp = sqrt(coresp(:,2).^2 + coresp(:,1).^2);
%     coresp = coresp(find((maskp<16.1).*(maskp>4.5)),:);
%     maskm = sqrt(coresm(:,2).^2 + coresm(:,1).^2); 
%     coresm = coresm(find((maskm<16.1).*(maskm>4.5)),:);
ip = (abs(coresp(:,3))<4); % & (abs(coresp(:,1).^2+coresp(:,2).^2)<100);
im = (abs(coresm(:,3))<4); % & (abs(coresm(:,1).^2+coresm(:,2).^2)<100);
    plot3(coresp(ip,1),coresp(ip,2),coresp(ip,3),'MarkerSize',16,'Marker','.','LineStyle','none','Clipping','off','Color',[0.5 0.5 0.5]);
    plot3(coresm(im,1),coresm(im,2),coresm(im,3),'MarkerSize',16,'Marker','.','LineStyle','none','Clipping','off','Color',[0.2 0.2 0.2]);
    xl=xlabel('X','FontSize',24);
lpos = get(xl,'Position');
lpos(1) = lpos(1) - 2;
set(xl,'Position',lpos);
    yl=ylabel('Y','FontSize',24);
lpos = get(yl,'Position');

lpos(2) = lpos(2) + 1;
set(yl,'Position',lpos);
    zl=zlabel('Z','FontSize',24);
        annotation(f,'textbox',...
        [0.596367521367525 0.781488457343255 0.217240808313939 0.0638361888494952],...
        'String',{sprintf('t = %0.3f s',(i*dt))},...
        'FontWeight','bold',...
        'FontSize',24,...
        'FitBoxToText','on',...
        'LineStyle','none',...
        'Color',[0 0 0]);
    
    % Create surf
patch(fv,'Parent',axes1,'LineStyle','none','FaceLighting','gouraud',...
    'FaceColor',[0.5 0.5 0.5],...
    'FaceAlpha',0.2,...
    'Clipping','off',...
    'EdgeAlpha',0.2);
zlpos = get(zl,'Position');
zlpos(1) = zlpos(1) + 1.2;
set(zl,'Position',zlpos);
% Create light
%camlight
light('Parent',axes1,...
    'Position',[-0.586743862956797 -0.0533583438964443 0.808012701892219]);
    set(f, 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 5]);
    print(f,'-dpng',sprintf('snapshots_m_1_s_2_om_0_25/corei_%05d.png',i), '-r150');
end