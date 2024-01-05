function ret = lz(obj,phi,r0)
% Calculates action of the Lz operator on the arbitrary state phi
% Optioanal argumant r0 allows to set the coordinates of the center

    if(nargin<=2)
        r0=[0,0];
    end
	ret = -1i*((obj.mesh.x-r0(1)).*obj.derivy(phi) - (obj.mesh.y-r0(2)).*obj.derivx(phi));
end
