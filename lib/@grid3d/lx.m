function ret = lx(obj,phi)
	ret = -1i*(obj.mesh.y.*obj.derivz(phi) - obj.mesh.z.*obj.derivy(phi));
end
