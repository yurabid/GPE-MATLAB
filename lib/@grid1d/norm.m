function n = norm( obj, psi )
%  NORM - Norm of wavefunction.
%
%  Usage for obj = grid3d :
%    n = norm( obj, psi )
%  Input
%    psi    :  wavefunction
%  Output
%    n      :  norm of wavefunction

n = sqrt(obj.integrate(abs(psi).^2));
