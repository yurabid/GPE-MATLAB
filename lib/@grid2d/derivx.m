function ret = derivx(obj,phi)
	ret = (-circshift(phi,[0 2])+8*circshift(phi,[0 1])-8*circshift(phi,[0 -1])+circshift(phi,[0 -2]))./(12*(obj.x(2)-obj.x(1)));
end
