function ret = derivy(obj,phi)
	ret = (-circshift(phi,[2 0])+8*circshift(phi,[1 0])-8*circshift(phi,[-1 0])+circshift(phi,[-2 0]))./(12*(obj.y(2)-obj.y(1)));
end
