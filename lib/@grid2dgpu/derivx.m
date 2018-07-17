function ret = derivx(obj,phi)
    nx = obj.nx;
    ind = [nx-1,nx,1:nx,1,2];
% 	ret = (-circshift(phi,[0 2])+8*circshift(phi,[0 1])-8*circshift(phi,[0 -1])+circshift(phi,[0 -2]))./(12*(obj.x(2)-obj.x(1)));
    ret = (-phi(:,ind(1:nx)) + 8*phi(:,ind(2:nx+1)) - 8*phi(:,ind(4:nx+3)) + phi(:,ind(5:nx+4)))/(12*(obj.x(2)-obj.x(1)));
end
