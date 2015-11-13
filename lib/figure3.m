function figure3(rx, ry, cdata, tarr, tind, Y1, Y2, maxx)

% Create figure
figure1 = gcf;
clf(figure1,'reset');
%colormap('jet');
fn = 'Arial';
[X,Y] = meshgrid(rx,ry);
[coresp, coresm] = detect_core(cdata,X,Y);
% Create axes
axes1 = axes('Parent',figure1,'YTick',[-10 0 10],...
    'Position',[0.085 0.528 0.31275 0.417],...
    'Layer','top', 'FontName',fn, 'LineWidth', 1);
 xlim(axes1,[rx(1) rx(end)]);
 ylim(axes1,[ry(1) ry(end)]);
box(axes1,'on');
hold(axes1,'all');
caxis(axes1,[0 maxx]);

% Create image
image(rx,ry,abs(cdata),'Parent',axes1,'CDataMapping','scaled');
%  if(tind < 785)
     plot(coresp(:,1),coresp(:,2),'Parent',axes1,'MarkerSize',4,'Marker','+','LineStyle','none','Color',[1 1 1], 'LineWidth', 1);
     plot(coresm(:,1),coresm(:,2),'Parent',axes1,'MarkerSize',4,'Marker','o','LineStyle','none','Color',[1 0 0], 'LineWidth', 1);
%  end
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
 xlim(axes2,[rx(1) rx(end)]);
 ylim(axes2,[ry(1) ry(end)]);
box(axes2,'on');
hold(axes2,'all');
caxis(axes2,[-pi pi]);

% Create image
image(rx, ry, angle(cdata),'Parent',axes2,'CDataMapping','scaled');
%  if(tind < 785)
 	plot(coresp(:,1),coresp(:,2),'Parent',axes2,'MarkerSize',4,'Marker','+','LineStyle','none','Color',[1 1 1], 'LineWidth', 1);
    plot(coresm(:,1),coresm(:,2),'Parent',axes2,'MarkerSize',4,'Marker','o','LineStyle','none','Color',[1 0 0], 'LineWidth', 1);
%  end
% Create xlabel
xlabel('x', 'FontAngle', 'italic', 'FontSize',13, 'FontName',fn);

% Create ylabel
yl=ylabel('y', 'FontAngle', 'italic', 'FontSize',13, 'FontName',fn);
ylpos = get(yl,'Position');
ylpos(1) = ylpos(1) + 2;
set(yl,'Position',ylpos);
title('Phase: Arg(\Psi(x,y,0))', 'FontSize',12, 'FontName',fn);


if(Y2>0)
	axes4 = axes('Parent',figure1,'Position',[0.085 0.085 0.835329341317365 0.247835497835498], 'FontName',fn, 'LineWidth', 1, 'XAxisLocation', 'top', 'XTick', []);
	xlim(axes4,[0 tarr(end)]);
	ylim(axes4,[0 Y2*2.9]);
	%line([0 0.5 1 1.5 tarr(end)], [0 Y2 Y2 0 Y2*(tarr(end)-1.5)/0.5],'Parent',axes4, 'LineWidth', 1, 'Color', [0 0.5 0]);
	line([0 0.5 1 1.5 tarr(end)], [0 Y2 Y2 0 0],'Parent',axes4, 'LineWidth', 1.5, 'Color', [0 0.5 0], 'LineStyle', '--');
	ylabel('U(t), kHz', 'FontAngle', 'italic', 'FontSize',13, 'FontName',fn, 'Color', [0 0.5 0]);
    axes3 = axes('Parent',figure1,'YGrid','on',...
        'Position',[0.085 0.085 0.835329341317365 0.247835497835498], 'FontName',fn, 'LineWidth', 1,'Color','none', 'YAxisLocation', 'right');
else
    axes3 = axes('Parent',figure1,'YGrid','on',...
     'Position',[0.085 0.085 0.835329341317365 0.247835497835498], 'FontName',fn, 'LineWidth', 1,'Color','none');    
end
% Create axes
% axes3 = axes('Parent',figure1,'YGrid','on',...
%     'Position',[0.085 0.085 0.835329341317365 0.247835497835498], 'FontName',fn, 'LineWidth', 1,'Color','none');
% box(axes3,'on');
hold(axes3,'all');
% if(max(Y1) > 5)
%     yl = ceil(max(Y1))+1;
% else
%     yl = ceil(max(Y1));
% end
% if(yl<10)
%     set(axes3,'YTick',[0 1 2 3 4 5 6 7 8 9 10]);
% elseif(yl<20)
%     set(axes3,'YTick',[0 2 4 6 8 10 12 14 16 18 20]);
% else
%     set(axes3,'YTick',[0 3 6 9 12 15 18 21 24 27 30]);
% end
%  ylim(axes3,[0 yl]);
xlim(axes3,[0 tarr(end)]);

% Create plot
line(tarr, Y1,'Parent',axes3, 'LineWidth', 1.5, 'Color', 'blue');

line(tarr(tind), Y1(tind),'Parent',axes3,'Marker','o','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','red');
mylim = ylim(axes3);
line([tarr(tind) tarr(tind)], mylim,'Parent',axes3, 'LineStyle', ':', 'LineWidth', 1,'Color','red');
% Create xlabel
xlabel('t, s', 'FontAngle', 'italic', 'FontSize',13, 'FontName',fn);

% Create ylabel
ylabel('L_z / N', 'FontAngle', 'italic', 'FontSize',13, 'FontName',fn, 'Color', 'blue');
%title('Angular momentum per particle / Height of the rotating barrier', 'FontSize',12, 'FontName',fn);
title('Angular momentum per particle', 'FontSize',12, 'FontName',fn);


% Create colorbar
colorbar(axes1,'Position',[0.41 0.528 0.016 0.417], 'FontName',fn, 'LineWidth', 1);

% Create colorbar
colorbar(axes2,'Position',[0.905 0.528 0.016 0.417],'YTick',[-3.14 0 3.14], 'FontName',fn, 'LineWidth', 1);

% Create textbox
annotation(figure1,'textbox',...
    [0.3 0.33 0.4 0.15],...
    'String',{sprintf('Time = %0.3f s', tarr(tind))},...
    'FontSize',16,...
    'HorizontalAlignment','center',...
    'FitBoxToText','off',...
    'LineStyle','none', 'FontName',fn);

% saveas(figure1,sprintf('combined_%05d.png',tind));
set(figure1, 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 6]);
print(figure1,'-dpng',sprintf('snapshots/combined_%05d.png',tind), '-r100');