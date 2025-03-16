function ret = lx(obj,phi,r0)
%% Calculates action of the Lx (angular momentum projection) operator on the arbitrary state phi
%% Optioanal argumant r0 allows to set the coordinates of the center

    if(nargin<=2)
        r0=[0,0,0];
    end
	ret = -1i*((obj.mesh.y-r0(2)).*obj.derivz(phi) - (obj.mesh.z-r0(3)).*obj.derivy(phi));
end
