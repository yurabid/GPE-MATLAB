function figure4(rx, ry, cdata, cdata2, tarr, tind, Y1, Y2, maxx)

% Create figure
figure1 = gcf;
clf(figure1,'reset');
fn = 'Arial';
[X,Y] = meshgrid(rx,ry);
[coresp, coresm] = detect_core(cdata,X,Y);
[coresp2, coresm2] = detect_core(cdata2,X,Y);
% Create axes
axes1 = axes('Parent',figure1,'YTick',[-10 0 10],...
    'Position',[0.085 0.528 0.31275 0.417],...
    'Layer','top', 'FontName',fn, 'LineWidth', 1);
 xlim(axes1,[-19.2 19.2]);
 ylim(axes1,[-19.2 19.2]);
box(axes1,'on');
hold(axes1,'all');
caxis(axes1,[0 maxx]);

% Create image
image(rx,ry,abs(cdata),'Parent',axes1,'CDataMapping','scaled');
	plot(coresp(:,1),coresp(:,2),'Parent',axes1,'MarkerSize',6,'Marker','+','LineStyle','none','Color',[1 1 1], 'LineWidth', 1.5);
    plot(coresm(:,1),coresm(:,2),'Parent',axes1,'MarkerSize',6,'Marker','o','LineStyle','none','Color',[1 0 0], 'LineWidth', 1.5);
% Create xlabel
xlabel('x', 'FontAngle', 'italic', 'FontSize',13, 'FontName',fn);

% Create ylabel
yl=ylabel('y', 'FontAngle', 'italic', 'FontSize',13, 'FontName',fn);
ylpos = get(yl,'Position');
ylpos(1) = ylpos(1) + 2;
set(yl,'Position',ylpos);
title('|\Psi(x,y,0)|', 'FontSize',12, 'FontName',fn);

% Create axes
axes2 = axes('Parent',figure1,'YTick',[-10 0 10],...
    'Position',[0.58 0.528 0.31275 0.417],...
    'Layer','top', 'FontName',fn, 'LineWidth', 1);
 xlim(axes2,[-19.2 19.2]);
 ylim(axes2,[-19.2 19.2]);
box(axes2,'on');
hold(axes2,'all');
%caxis(axes2,[-pi pi]);
caxis(axes2,[-pi pi]);

% Create image
image(rx, ry, angle(cdata),'Parent',axes2,'CDataMapping','scaled');
	plot(coresp(:,1),coresp(:,2),'Parent',axes2,'MarkerSize',6,'Marker','+','LineStyle','none','Color',[1 1 1], 'LineWidth', 1.5);
    plot(coresm(:,1),coresm(:,2),'Parent',axes2,'MarkerSize',6,'Marker','o','LineStyle','none','Color',[1 0 0], 'LineWidth', 1.5);
% Create xlabel
xlabel('x', 'FontAngle', 'italic', 'FontSize',13, 'FontName',fn);

% Create ylabel
yl=ylabel('y', 'FontAngle', 'italic', 'FontSize',13, 'FontName',fn);
ylpos = get(yl,'Position');
ylpos(1) = ylpos(1) + 2;
set(yl,'Position',ylpos);
title('Phase: Arg(\Psi(x,y,0))', 'FontSize',12, 'FontName',fn);


% Create axes
axes3 = axes('Parent',figure1,'YTick',[-10 0 10],...
    'Position',[0.085 0.085 0.31275 0.417],...
    'Layer','top', 'FontName',fn, 'LineWidth', 1);
 xlim(axes3,[-19.2 19.2]);
 ylim(axes3,[-19.2 19.2]);
box(axes3,'on');
hold(axes3,'all');
caxis(axes3,[0 maxx]);

% Create image
image(rx,ry,abs(cdata2),'Parent',axes3,'CDataMapping','scaled');
	plot(coresp2(:,1),coresp2(:,2),'Parent',axes3,'MarkerSize',6,'Marker','+','LineStyle','none','Color',[1 1 1], 'LineWidth', 1.5);
    plot(coresm2(:,1),coresm2(:,2),'Parent',axes3,'MarkerSize',6,'Marker','o','LineStyle','none','Color',[1 0 0], 'LineWidth', 1.5);
% Create xlabel
xlabel('y', 'FontAngle', 'italic', 'FontSize',13, 'FontName',fn);

% Create ylabel
yl=ylabel('z', 'FontAngle', 'italic', 'FontSize',13, 'FontName',fn);
ylpos = get(yl,'Position');
ylpos(1) = ylpos(1) + 2;
set(yl,'Position',ylpos);
%title('|\Psi(x,y,0)|', 'FontSize',12, 'FontName',fn);

% Create axes
axes4 = axes('Parent',figure1,'YTick',[-10 0 10],...
    'Position',[0.58 0.085 0.31275 0.417],...
    'Layer','top', 'FontName',fn, 'LineWidth', 1);
 xlim(axes4,[-19.2 19.2]);
 ylim(axes4,[-19.2 19.2]);
box(axes4,'on');
hold(axes4,'all');
caxis(axes4,[-pi pi]);

% Create image
image(rx, ry, angle(cdata2),'Parent',axes4,'CDataMapping','scaled');
	plot(coresp2(:,1),coresp2(:,2),'Parent',axes4,'MarkerSize',6,'Marker','+','LineStyle','none','Color',[1 1 1], 'LineWidth', 1.5);
    plot(coresm2(:,1),coresm2(:,2),'Parent',axes4,'MarkerSize',6,'Marker','o','LineStyle','none','Color',[1 0 0], 'LineWidth', 1.5);
% Create xlabel
xlabel('y', 'FontAngle', 'italic', 'FontSize',13, 'FontName',fn);

% Create ylabel
yl=ylabel('z', 'FontAngle', 'italic', 'FontSize',13, 'FontName',fn);
ylpos = get(yl,'Position');
ylpos(1) = ylpos(1) + 2;
set(yl,'Position',ylpos);
%title('Phase: Arg(\Psi(x,y,0))', 'FontSize',12, 'FontName',fn);






% Create colorbar
colorbar(axes1,'Position',[0.41 0.528 0.016 0.417], 'FontName',fn, 'LineWidth', 1);
colorbar(axes2,'Position',[0.905 0.528 0.016 0.417],'YTick',[-3.14 0 3.14], 'FontName',fn, 'LineWidth', 1);

% Create colorbar
colorbar(axes3,'Position',[0.41 0.085 0.016 0.417], 'FontName',fn, 'LineWidth', 1);
colorbar(axes4,'Position',[0.905 0.085 0.016 0.417],'YTick',[-3.14 0 3.14], 'FontName',fn, 'LineWidth', 1);

% Create textbox
annotation(figure1,'textbox',...
    [0.3 0.33 0.4 0.15],...
    'String',{sprintf('t = %0.3f s', tarr(tind))},...
    'FontSize',16,...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'LineStyle','none', 'FontName',fn);

% saveas(figure1,sprintf('combined_%05d.png',tind));
set(figure1, 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 6]);
print(figure1,'-dpng',sprintf('snapshots/combined_%05d.png',tind), '-r150');
