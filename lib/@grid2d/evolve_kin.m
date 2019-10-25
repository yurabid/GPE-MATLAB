function phi = evolve_kin(obj,phi,dt,is_inverse,omega)
    if(nargin<4)
        is_inverse = 0;
    end
    if(nargin<5)
        omega=0;
    end
    if(is_inverse == 0)
        ekk = exp(-(obj.kx.^2/2+omega*obj.mesh.y)*dt);
        phi = ifft(ekk.*fft(phi,[],1),[],1);
    end
    ekk = exp(-(obj.ky.^2/2-omega*obj.mesh.x)*dt);
    phi = ifft(ekk.*fft(phi,[],2),[],2);
    if(is_inverse ~= 0)
        ekk = exp(-(obj.kx.^2/2+omega*obj.mesh.y)*dt);
        phi = ifft(ekk.*fft(phi,[],1),[],1);
    end
end