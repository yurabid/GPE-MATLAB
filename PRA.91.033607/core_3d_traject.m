f = figure;
dt = 4*1.29366e-3;
% VG = gather(V(:,:,16));
for i = 1:1000
    load(sprintf('core_%05d.mat',i));
    clf(f,'reset');
    axes1 = axes('Parent',f);
    xlim(axes1,[-13 13]);
    ylim(axes1,[-13 13]);
    zlim(axes1,[-8 8]);
    view(axes1,[-40.5 24]);
    box(axes1,'on');
    grid(axes1,'on');
    hold(axes1,'all');
    %plot3(cores(:,1),cores(:,2),cores(:,3),'MarkerSize',2,'Marker','o','LineStyle','none');
    maskp = sqrt(coresp(:,2).^2 + coresp(:,1).^2);
    coresp = coresp(find((maskp<16.1).*(maskp>4.5)),:);
    maskm = sqrt(coresm(:,2).^2 + coresm(:,1).^2); 
    coresm = coresm(find((maskm<16.1).*(maskm>4.5)),:);
    plot3(coresp(:,1),coresp(:,2),coresp(:,3),'MarkerSize',6,'Marker','.','LineStyle','none','Color',[1 0 0]);
    plot3(coresm(:,1),coresm(:,2),coresm(:,3),'MarkerSize',6,'Marker','.','LineStyle','none','Color',[0 0 1]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
        annotation(f,'textbox',...
        [0.7 0.8 0.21092023633678 0.127424242424244],...
        'String',{sprintf('t = %0.2f',(i*dt))},...
        'FontWeight','bold',...
        'FontSize',20,...
        'FitBoxToText','on',...
        'LineStyle','none',...
        'Color',[0 0 0]);
    
    % Create surf
patch(fv2,'Parent',axes1,'LineStyle','none',...
    'FaceColor',[0.5 0.5 0.5],...
    'FaceAlpha',0.3,...
    'EdgeAlpha',0.3);

% Create light
%camlight
light('Parent',axes1,...
    'Position',[-0.586743862956797 -0.0533583438964443 0.808012701892219]);
    
    saveas(f,sprintf('core_%05d.png',i));
end