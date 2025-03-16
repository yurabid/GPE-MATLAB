function ret = lap(obj,phi,omega)
% Applies the kinetic energy operator to the arbitrary state phi
% Optional argument omega allows to do the calculation in a rotating frame

if(nargin <= 2 || omega==0)
	ret = obj.ifft(obj.kk.*obj.fft(phi));
else
	ret = obj.ifftx((obj.kx.^2/2-obj.kx.*obj.mesh.y*omega).*obj.fftx(phi)) + ...
		obj.iffty((obj.ky.^2/2+obj.ky.*obj.mesh.x*omega).*obj.ffty(phi)) + ...
                obj.ifftz((obj.kz.^2/2).*obj.fftz(phi));
end
end
