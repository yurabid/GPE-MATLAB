function ret = lap(obj,phi,omega)
if(nargin <= 2)
	ret = obj.ifft(obj.kk.*obj.fft(phi));
else
	ret = obj.ifftx((obj.kx.^2/2-obj.kx.*obj.mesh.y*omega).*obj.fftx(phi)) + ...
		obj.iffty((obj.ky.^2/2+obj.ky.*obj.mesh.x*omega).*obj.ffty(phi));
end
end
