function ret = derivx(obj,phi)
	ret = (-circshift(phi,2)+8*circshift(phi,1)-8*circshift(phi,-1)+circshift(phi,-2))./(12*(obj.x(2)-obj.x(1)));
end
