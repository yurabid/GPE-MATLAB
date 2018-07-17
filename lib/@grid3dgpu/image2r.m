function [a1, a2] = image2r(obj,phi1,phi2,varargin)
% image2 - Make a two-panel image of density and phase.
%
%  Usage :
%    [a1, a2] = task.image2(data)
%    [a1, a2] = task.image2(data,figure)
%    [a1, a2] = task.image2(data,figure,bbox)
%    [a1, a2] = task.image2(data,figure,bbox,....) - pass additional
%    options to axes and colorbar functions
%    task.image2() - reset color scale for density plot
%  Input
%    data    :  2D or 3D data array
%    figure  :  figure handle, default - current figure
%    bbox    :  bounding box of added panels, default [0 0 1 1]
%  Output
%    a1, a2   :  axes handles of two created panels

persistent maxx

if(nargin==1 || numel(phi1) <= 1)
    clear maxx;
    return;
end
opts = {};
if (nargin>3)
    f = varargin{1};
else
    f = gcf;
    clf(f,'reset');
    set(f,'Position',[100 200 700 280]);
end
if (nargin>4)
    bbox = varargin{2};
else
    bbox = [0 0 1 1];
end
if (nargin>5)
    opts = varargin(3:end);
end
bbwidth = bbox(3)-bbox(1);
bbheight = bbox(4)-bbox(2);

if (ndims(phi1) == 3)
    slice1 = phi1(:,:,obj.nz/2);
    slice2 = phi2(:,:,obj.nz/2);
else
    slice1 = phi1;
    slice2 = phi2;
end
cdata1 = gather(slice1);
cdata2 = gather(slice2);
rx = gather(obj.x);
ry = gather(obj.y);

if (isempty(maxx))
    maxx = max(cdata1(:));
end

a1 = axes('Parent',f,'Position',[bbox(1)+bbwidth*0.1 bbox(2)+bbheight*0.15 bbwidth*0.32 bbheight*0.72], 'LineWidth', 1, opts{:});
a2 = axes('Parent',f,'Position',[bbox(1)+bbwidth*0.59 bbox(2)+bbheight*0.15 bbwidth*0.32 bbheight*0.72], 'LineWidth', 1, opts{:});

box(a1,'on');
hold(a1,'all');
caxis(a1,[0 maxx]);
xlim(a1,[rx(1) rx(end)]);
ylim(a1,[ry(1) ry(end)]);

box(a2,'on');
hold(a2,'all');
caxis(a2,[0 maxx]);
xlim(a2,[rx(1) rx(end)]);
ylim(a2,[ry(1) ry(end)]);

image(rx,ry,cdata1,'Parent',a1,'CDataMapping','scaled');
image(rx,ry,cdata2,'Parent',a2,'CDataMapping','scaled');

%colorbar(a1,'Position',[bbox(1)+bbwidth*0.435 bbox(2)+bbheight*0.15 bbwidth*0.02 bbheight*0.72], 'LineWidth', 1, opts{:});
%colorbar(a2,'Position',[bbox(1)+bbwidth*0.925 bbox(2)+bbheight*0.15 bbwidth*0.02 bbheight*0.72], 'LineWidth', 1, opts{:});

end
