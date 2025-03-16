function ret = derivz(obj,phi)
    ret = obj.deriv1(phi,3,obj.dx);
	% ret = (-circshift(phi,[0 0 2])+8*circshift(phi,[0 0 1])-8*circshift(phi,[0 0 -1])+circshift(phi,[0 0 -2]))./(12*(obj.z(2)-obj.z(1)));
end
