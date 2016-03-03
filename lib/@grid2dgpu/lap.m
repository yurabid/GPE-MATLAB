function ret = lz(obj,phi)
	ret = obj.ifft(obj.kk.*obj.fft(phi));
end
