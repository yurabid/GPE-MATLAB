function int = integrate( obj, fun )
%  INTEGRATE - Integrate function over grid domain.
%
%  Usage for obj = grid3d :
%    int = integrate( obj, fun )
%  Input
%    fun    :  integrand
%  Output
%    int    :  integrated function
f = fun.*obj.mesh.r;
int = sum(f(:)) * obj.weight;