function I = sommerfelds(kx,x,k1,k2,eps1,eps2)
% (kx,x,k1,k2,eps1,eps2) returns integrand of a Sommerfield integral,
% used for calculations of Nevels-Michalski's magnetic field equations.
%
% kx is the dummy variable used by the integration routine
% See PQintegral for further information.
y_int = imag(k2^2)/(2*k1);
kz1 = sqrt( k1.*k1 - (kx.^ 2));
kz2 = sqrt( k2.*k2 - (kx.^ 2));
D = kz1./eps1 + kz2./eps2;
G = 1./D;
I = G .* exp(-1i*kx*x);

end