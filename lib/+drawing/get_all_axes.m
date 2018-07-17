function axs=get_all_axes(fig)

axs = [];

if(nargin==0)
    fig = gcf;
end

c = get(fig,'Children');

for i=1:length(c)
    if(isa(c(i),'matlab.graphics.axis.Axes'))
        axs = [c(i), axs];
    end
end
