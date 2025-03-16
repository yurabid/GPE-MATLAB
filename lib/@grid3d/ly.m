function ret = ly(obj,phi,r0)
%% Calculates action of the Ly (angular momentum projection) operator on the arbitrary state phi
%% Optioanal argumant r0 allows to set the coordinates of the center

    if(nargin<=2)
        r0=[0,0,0];
    end
	ret = -1i*((obj.mesh.z-r0(3)).*obj.derivx(phi) - (obj.mesh.x-r0(1)).*obj.derivz(phi));
end
