function ret = derivz(obj,phi)
	ret = (-circshift(phi,[0 0 2])+8*circshift(phi,[0 0 1])-8*circshift(phi,[0 0 -1])+circshift(phi,[0 0 -2]))./(12*(obj.z(2)-obj.z(1)));
end
