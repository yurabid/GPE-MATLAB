function [s,v] = velocity_field3d(phi,r,cphase)
% Calculate the velocity field on a 3D wave function
% IN:
%   - phi: 3D wave function
%   - r: grid of coordinate ponts (same in each direction)
%   - cphase: the angle on a plane at which the stramlines will begin
% OUT:
%   - v: scalar velocity field (absolute values)
%   - s: streamlines of the vector velocity field
	h = r(2)-r(1);
	[sx,sy,sz] = gradient(phi);
	[X,Y,Z] = meshgrid(r,r,r);
	
	vx = -imag(phi.*conj(sx) - conj(phi).*sx)./(abs(phi).^2).*(abs(phi)>1e-4)./(2*h);
	vy = -imag(phi.*conj(sy) - conj(phi).*sy)./(abs(phi).^2).*(abs(phi)>1e-4)./(2*h);
	vz = -imag(phi.*conj(sz) - conj(phi).*sz)./(abs(phi).^2).*(abs(phi)>1e-4)./(2*h);
	
	v = sqrt(vx.^2+vy.^2+vz.^2);
	%startr = [ -6:-0.7:-16 6:0.7:16];
	startr = [4];
	starty = startr*sin(cphase);
	startx = startr*cos(cphase);
	startz = [0];
	s = stream3(X,Y,Z,vx,vy,vz,startx,starty,startz,[0.2, 3000]);
	%s2 = stream3(X,Y,Z,-vx,-vy,-vz,startx,starty,startz,[0.2, 2000]);
	%s = [s1 s2];
