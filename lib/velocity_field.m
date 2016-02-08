function [s,v, varargout] = velocity_field(slice,r,cphase,slice0)
% Calculate the velocity field on a 2D slice of the wave function
% IN:
%   - slice: 2D slice of the wave function
%   - r: grid of coordinate ponts (same in each direction)
%   - cphase: the angle on a plane at which the stramlines will begin
% OUT:
%   - v: scalar velocity field (absolute values)
%   - s: streamlines of the vector velocity field
	h = r(2)-r(1);
	[sx,sy] = gradient(slice);
	[X,Y] = meshgrid(r,r);
	%vx0 = -imag(slice0.*conj(sx) - conj(slice0).*sx)./(abs(slice0).^2).*(abs(slice0)>1e-4)./(2*h);
	%vy0 = -imag(slice0.*conj(sy) - conj(slice0).*sy)./(abs(slice0).^2).*(abs(slice0)>1e-4)./(2*h);
	
	vx = -imag(slice.*conj(sx) - conj(slice).*sx)./(abs(slice).^2).*(abs(slice)>1e-6)./(2*h);% - vx0;
	vy = -imag(slice.*conj(sy) - conj(slice).*sy)./(abs(slice).^2).*(abs(slice)>1e-6)./(2*h);% - vy0;
	
	v = sqrt(vx.^2+vy.^2);
	startr = [  6:0.7:16];
	starty = [startr*sin(cphase) -startr*sin(cphase)];
	startx = [startr*cos(cphase) startr*cos(cphase)];
	s1 = stream2(X,Y,vx,vy,startx,starty,[0.2, 1000]);
	s2 = stream2(X,Y,-vx,-vy,startx,starty,[0.2, 1000]);
	s = [s1 s2];
    if (nargout == 4)
        varargout{1} = vx;
        varargout{2} = vy;
    end
