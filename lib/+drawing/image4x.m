function [axs,cbars]=image4x(x, y, cdata, cdata2, cdata3, cdata4, showcb, hideshow, forcereset, handle)

% Create figure
if(nargin>9)
    figure1 = handle;
else
    figure1 = figure(4000);
end
if(nargin<8 || hideshow>0)
    set(figure1, 'Visible', 'off');
end
rx = gather(x);
ry = gather(y);
cdata = gather(real(cdata)); cdata2=gather(real(cdata2)); cdata3=gather(real(cdata3)); cdata4=gather(real(cdata4));
axs = [];
cbars = [];
if(nargin<9)
    forcereset = 0;
end
if(forcereset==0)
    axs = drawing.get_all_axes(figure1);
    if(length(axs)~=4)
        forcereset=1;
    end
end
if(forcereset==0)
%     axs = drawing.get_all_axes(figure1);
    cla(axs(1));cla(axs(2));cla(axs(3));cla(axs(4));
else

clf(figure1,'reset');
fn = 'Arial';
set(figure1, 'DefaultTextInterpreter', 'latex');
% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.085 0.58 0.35 0.35],...
    'Layer','top', 'LineWidth', 1);
 xlim(axes1,[rx(1) rx(end)]);
 ylim(axes1,[ry(1) ry(end)]);
box(axes1,'on');
xlabel('$x$', 'FontSize',16);
yl=ylabel('$y$', 'FontSize',16);
ylpos = get(yl,'Position');
ylpos(1) = ylpos(1) + 2;
set(yl,'Position',ylpos);
% title('$|\Psi(x,y,0)|$', 'FontSize',12, 'FontName',fn);

% Create axes
axes2 = axes('Parent',figure1,...
    'Position',[0.58 0.58 0.35 0.35],...
    'Layer','top', 'LineWidth', 1);
 xlim(axes2,[rx(1) rx(end)]);
 ylim(axes2,[ry(1) ry(end)]);
box(axes2,'on');
xlabel('$x$', 'FontSize',16);
yl=ylabel('$y$', 'FontSize',16);
ylpos = get(yl,'Position');
ylpos(1) = ylpos(1) + 2;
set(yl,'Position',ylpos);
% title('Phase: Arg(\Psi(x,y,0))', 'FontSize',12, 'FontName',fn);


% Create axes
axes3 = axes('Parent',figure1,...
    'Position',[0.085 0.085 0.35 0.35],...
    'Layer','top', 'LineWidth', 1);
 xlim(axes3,[rx(1) rx(end)]);
 ylim(axes3,[ry(1) ry(end)]);
box(axes3,'on');
xlabel('$x$', 'FontSize',16);
yl=ylabel('$y$', 'FontSize',16);
ylpos = get(yl,'Position');
ylpos(1) = ylpos(1) + 2;
set(yl,'Position',ylpos);
%title('|\Psi(x,y,0)|', 'FontSize',12, 'FontName',fn);

% Create axes
axes4 = axes('Parent',figure1,...
    'Position',[0.58 0.085 0.35 0.35],...
    'Layer','top', 'LineWidth', 1);
 xlim(axes4,[rx(1) rx(end)]);
 ylim(axes4,[ry(1) ry(end)]);
box(axes4,'on');
xlabel('$x$', 'FontSize',16);
yl=ylabel('$y$',  'FontSize',16);
ylpos = get(yl,'Position');
ylpos(1) = ylpos(1) + 2;
set(yl,'Position',ylpos);
%title('Phase: Arg(\Psi(x,y,0))', 'FontSize',12, 'FontName',fn);

set(axes1, 'XTickMode', 'manual');
set(axes1, 'YTickMode', 'manual');
set(axes1, 'XTickLabelMode', 'manual');
set(axes1, 'YTickLabelMode', 'manual');
set(axes2, 'XTickMode', 'manual');
set(axes2, 'YTickMode', 'manual');
set(axes2, 'XTickLabelMode', 'manual');
set(axes2, 'YTickLabelMode', 'manual');
set(axes3, 'XTickMode', 'manual');
set(axes3, 'YTickMode', 'manual');
set(axes3, 'XTickLabelMode', 'manual');
set(axes3, 'YTickLabelMode', 'manual');
set(axes4, 'XTickMode', 'manual');
set(axes4, 'YTickMode', 'manual');
set(axes4, 'XTickLabelMode', 'manual');
set(axes4, 'YTickLabelMode', 'manual');
set(axes1, 'CLimMode', 'manual');
set(axes2, 'CLimMode', 'manual');
set(axes3, 'CLimMode', 'manual');
set(axes4, 'CLimMode', 'manual');
set(axes1, 'CLim', [min(cdata(:)), max(real(cdata(:)))]);
set(axes2, 'CLim', [min(cdata2(:)), max(cdata2(:))]);
set(axes3, 'CLim', [min(cdata3(:)), max(cdata3(:))]);
set(axes4, 'CLim', [min(cdata4(:)), max(cdata4(:))]);
    
if(nargin>6 && showcb>0)
    % Create colorbar
    cbar1 = colorbar(axes1,'Position',[0.45 0.58 0.02 0.35], 'FontName',fn, 'LineWidth', 1);
    cbar2 = colorbar(axes2,'Position',[0.945 0.58 0.02 0.35],'FontName',fn, 'LineWidth', 1);
    cbar3 = colorbar(axes3,'Position',[0.45 0.085 0.02 0.35], 'FontName',fn, 'LineWidth', 1);
    cbar4 = colorbar(axes4,'Position',[0.945 0.085 0.02 0.35],'FontName',fn, 'LineWidth', 1);
    cbars = [cbar1,cbar2,cbar3,cbar4];
end
% saveas(figure1,sprintf('combined_%05d.png',tind));
cpos = get(figure1, 'pos');
cpos(3) = 700;
cpos(4) = 700;
set(figure1, 'pos', cpos);
axs = [axes1,axes2,axes3,axes4];
hold(axs(1),'on');hold(axs(2),'on');hold(axs(3),'on');hold(axs(4),'on');

end
% Create image
surface(rx,ry,zeros(size(cdata)),cdata,'Parent',axs(1),'CDataMapping','scaled','EdgeColor','none');
surface(rx,ry,zeros(size(cdata2)),cdata2,'Parent',axs(2),'CDataMapping','scaled','EdgeColor','none');
surface(rx,ry,zeros(size(cdata3)),cdata3,'Parent',axs(3),'CDataMapping','scaled','EdgeColor','none');
surface(rx,ry,zeros(size(cdata4)),cdata4,'Parent',axs(4),'CDataMapping','scaled','EdgeColor','none');

% set(figure1, 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 10]);
%print(figure1,'-dpng',sprintf('snapshots/combined_%05d.png',tind), '-r150');
if(nargin<8 || hideshow>0)
    set(figure1, 'Visible', 'on');
end