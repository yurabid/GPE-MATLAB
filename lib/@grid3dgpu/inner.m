function y = inner( obj, lhs, rhs )
%  INNER - Inner product <lhs,rhs> of functions lhs and rhs.
%
%  Usage for obj = grid3d :
%    y = inner( obj, lhs, rhs )
%  Input
%    lhs    :  left  hand side function
%    rhs    :  right hand side function
%  Output
%    y      :  integrated inner product

y = integrate( obj, conj( lhs ) .* rhs );
