f = figure;

axes1 = axes('Parent',f);
    xlim(axes1,[-13 13]);
    ylim(axes1,[-13 13]);
    zlim(axes1,[-10 10]);
    view(axes1,[-40.5 24]);
    box(axes1,'off');
    grid(axes1,'off');
    hold(axes1,'all');
    [X,Y,Z] = meshgrid(r,r, [-0.2 0.2]);
    fv3 = isosurface(X,Y,Z,repmat(VVV,1,1,2),12); 
    
% Create surf
patch(fv,'Parent',axes1,'LineStyle','none',...
    'FaceColor',[0.5 0.5 0.5],...
    'FaceAlpha',0.3,...
    'EdgeAlpha',0.3);

% Create surf
patch(fv3,'Parent',axes1,'LineStyle','none',...
    'FaceColor',[0 1 0],...
    'FaceAlpha',0.3,...
    'EdgeAlpha',0.3);

% Create light
%camlight
light('Parent',axes1,...
    'Position',[-0.586743862956797 -0.0533583438964443 0.808012701892219]);