function ret = lz(obj,phi)
	ret = -1i*(obj.mesh.x.*obj.derivy(phi) - obj.mesh.y.*obj.derivx(phi));
end
