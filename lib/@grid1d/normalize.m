function [ psi, n ] = normalize( obj, psi )
%  NORMALIZE - Normalize wavefunction.
%
%  Usage for obj = grid3d :
%    [ psi, n ] = norm( obj, psi )
%  Input
%    psi    :  wavefunction
%  Output
%    psi    :  normalized wavefunction
%    n      :  norm of wavefunction

%  norm
n = norm( obj, psi );
%  normalize wavefunction
psi = psi * diag( 1 ./ n );

