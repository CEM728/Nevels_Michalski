function a = somint(k1,k2,eps1,eps2)
% PQintegral(D,T,Z,mode) performs a Sommerfeld. integral
% mode=l: returns P(D,T,Z) for up-link direction
% mode=2: returns P(D,T,Z) for down-link direction
% mode=3: returns Q(D,T,Z) for either direction
% Z is normalised height of observation point
% D is normalised offset of observation point
% T is normalised depth of transmitter
TOL = 1e-8; % the default is IE-6
TRACE = 0;
XLIMIT = k1;
% a = quadl(@sommerfelds,k1 - 100*k1*1i,XLIMIT,TOL,TRACE,k1,k2,eps1,eps2);
a = quadl(@sommerfelds,k1 - 100*k1*1i, k1);
% Numerically evaluate integral; adaptive Lobatto quadrature.