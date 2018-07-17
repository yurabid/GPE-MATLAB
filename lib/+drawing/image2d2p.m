function axs=image2d2p(x, y, cdata, cdata2, showcores, showcb)

% Create figure
figure1 = figure(2020);
if(nargin<6)
   showcb=0; 
end
if(nargin<5)
    showcores=0;
end

set(figure1, 'Visible', 'off');
[axs,cbars]=drawing.image4x(x, y, abs(cdata), angle(cdata), abs(cdata2), angle(cdata2), showcb, 0, 0, figure1);

if(showcores>0)
    [X,Y] = meshgrid(x,y);
    [coresp, coresm] = detect_core(cdata,X,Y);
    [coresp2, coresm2] = detect_core(cdata2,X,Y);
	plot(coresp(:,1),coresp(:,2),'Parent',axs(1),'MarkerSize',6,'Marker','+','LineStyle','none','Color',[1 1 1], 'LineWidth', 1.5);
    plot(coresm(:,1),coresm(:,2),'Parent',axs(1),'MarkerSize',6,'Marker','o','LineStyle','none','Color',[1 1 1], 'LineWidth', 1.5);
	plot(coresp(:,1),coresp(:,2),'Parent',axs(2),'MarkerSize',6,'Marker','+','LineStyle','none','Color',[1 1 1], 'LineWidth', 1.5);
    plot(coresm(:,1),coresm(:,2),'Parent',axs(2),'MarkerSize',6,'Marker','o','LineStyle','none','Color',[1 1 1], 'LineWidth', 1.5);
	plot(coresp2(:,1),coresp2(:,2),'Parent',axs(3),'MarkerSize',6,'Marker','+','LineStyle','none','Color',[1 1 1], 'LineWidth', 1.5);
    plot(coresm2(:,1),coresm2(:,2),'Parent',axs(3),'MarkerSize',6,'Marker','o','LineStyle','none','Color',[1 1 1], 'LineWidth', 1.5);
	plot(coresp2(:,1),coresp2(:,2),'Parent',axs(4),'MarkerSize',6,'Marker','+','LineStyle','none','Color',[1 1 1], 'LineWidth', 1.5);
    plot(coresm2(:,1),coresm2(:,2),'Parent',axs(4),'MarkerSize',6,'Marker','o','LineStyle','none','Color',[1 1 1], 'LineWidth', 1.5);
end

if(showcb>0 && ~isempty(cbars))
   cb = cbars(2);
   set(cb,'Ticks',[-3.14 0 3.14]);
   set(cb,'TickLabels',{'-\pi','0','\pi'});
   cb = cbars(4);
   set(cb,'Ticks',[-3.14 0 3.14]);
   set(cb,'TickLabels',{'-\pi','0','\pi'});
end

set(figure1, 'Visible', 'on');
