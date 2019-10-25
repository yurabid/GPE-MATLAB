function axs=image1d1p(x, y, cdata, showcores, showcb)

% Create figure
figure1 = figure(1010);
if(nargin<6)
   showcb=0; 
end
if(nargin<5)
    showcores=0;
end

set(figure1, 'Visible', 'off');
[axs,cbars]=drawing.image2x(x, y, (abs(cdata).^2), angle(cdata), showcb, 0, 0, figure1);

if(showcores>0)
    [X,Y] = meshgrid(x,y);
    [coresp, coresm] = detect_core(cdata,X,Y);
	plot(coresp(:,1),coresp(:,2),'Parent',axs(1),'MarkerSize',6,'Marker','+','LineStyle','none','Color',[1 1 1], 'LineWidth', 1.5);
    plot(coresm(:,1),coresm(:,2),'Parent',axs(1),'MarkerSize',6,'Marker','o','LineStyle','none','Color',[1 1 1], 'LineWidth', 1.5);
	plot(coresp(:,1),coresp(:,2),'Parent',axs(2),'MarkerSize',6,'Marker','+','LineStyle','none','Color',[1 1 1], 'LineWidth', 1.5);
    plot(coresm(:,1),coresm(:,2),'Parent',axs(2),'MarkerSize',6,'Marker','o','LineStyle','none','Color',[1 1 1], 'LineWidth', 1.5);
end

if(showcb>0 && ~isempty(cbars))
   cb = cbars(2);
   set(cb,'Ticks',[-3.14 0 3.14]);
   set(cb,'TickLabels',{'-\pi','0','\pi'});
end

set(figure1, 'Visible', 'on');
