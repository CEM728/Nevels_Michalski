function a = sommerfeld(x,D,T,Z,mode)
% (x,D,T,Z,mode) returns integrand of a Sommerfield integral,
% used for calculations of Wait's magnetic field equations.
%
% x is the dummy variable used by the integration routine
% See PQintegral for further information.
u = sqrt( x.*x + (T.^ 2)*1i );
if T == 0 % special case
    a = (x) .* exp(-u) .* exp(-x.*Z) /2;
else
    a = (x.^2) .* exp(-u) .* exp(-x.*Z) ./(x + u);
end
if mode == 1 % uplink P field
    a = x .* besselj(1, x .*D) .*a;
elseif mode == 2
    a = u .* besselj(1, x .*D) .*a;
elseif mode == 3
    a = x .* besselj(0, x .*D) .*a;
else
    a = O;
end