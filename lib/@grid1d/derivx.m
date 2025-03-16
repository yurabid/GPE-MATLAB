function ret = derivx(obj,phi)
    ret = obj.deriv1(phi,1,obj.dx);
	% ret = (-circshift(phi,2)+8*circshift(phi,1)-8*circshift(phi,-1)+circshift(phi,-2))./(12*(obj.x(2)-obj.x(1)));
end
