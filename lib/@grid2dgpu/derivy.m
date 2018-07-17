function ret = derivy(obj,phi)
    ny = obj.ny;
    ind = [ny-1,ny,1:ny,1,2];
% 	ret = (-circshift(phi,[2 0])+8*circshift(phi,[1 0])-8*circshift(phi,[-1 0])+circshift(phi,[-2 0]))./(12*(obj.y(2)-obj.y(1)));
    ret = (-phi(ind(1:ny),:) + 8*phi(ind(2:ny+1),:) - 8*phi(ind(4:ny+3),:) + phi(ind(5:ny+4),:))/(12*(obj.y(2)-obj.y(1)));
end
