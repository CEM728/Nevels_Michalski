clear all
close all

lambda = 633e-9; % Green light wavelength
eps_silver =  -18.295 - 1i*0.48085; % Johnson & Christy,1972 (refractiveindex.info) at 633 nm
load em_constants.mat % Contains varepsilon, mu and c
omega = 2*pi*c/lambda; % angular frequency

k_air = 2*pi/lambda; % propagation constant of air
% k_air = omega*sqrt(
k_silver = omega* sqrt(mu_0*(epsilon_0*eps_silver)); % propagation constant of silver
kxx = linspace(-1e4*k_air,1e4*k_air,1e4);
kxy = linspace(-1e4*k_air,1e4*k_air,1e4);
kx = meshgrid(kxx,kxy); % [real(kx),imag(kx)]
dif = abs(kxx - k_air);
pole_location = find(dif == min(dif)); % Index in kxx with the nearest value of k_air

x = linspace(1e-2*lambda,1e4*lambda,10); % Space vector definition
H = zeros(length(x)); % Initialize the magnetic field vector
%
%% Contour Definition
%
% kx_real = k_air;   kx_im_start = -k_air*100;   eps10 = .1; % Upward path
% s = chebfun('s',[0 1]);                 % dummy variable
% c = [ (k_air - eps10) + kx_im_start*1i ...
%     (k_air - eps10) + eps*1i ...
%     (k_air + eps10) + eps*1i ...
%     (k_air + eps10) + kx_im_start*1i ]; % Path points
% r = sqrt(kx_im_start^2 + eps10^2); % find the radius of the large enlosing circle
% alpha = atan(eps10/r); % arc till the starting point
% c0 = chebfun('9.926043139304243e+06 + eps*1e10*exp(-1i*s)',[-pi 0]);  % Semicircle above the branch-point
% c1 = chebfun('9.926043139304243e+06 + (9.926043139304243e+08)*exp(1i*s)',[3*pi/2+alpha (3*pi/2+2*pi-alpha)]);  % big cirle to infinity
%
%% Make a closed contour 
%
% kx = join( c(1) + s*(c(2)-c(1)), ...     % left branch going up
%     c0, ...                        % semi-circle surrounding the branch point
%     c(3) + s*(c(4)-c(3)), ...     % right branch going down
%     c1 );
%
%% Integrate point by point
%
% for i = 1 : length (x)
%     G = @(kx) 1./((sqrt(k_air^2 - kx.^2)./epsilon_0 + sqrt(k_silver^2 - kx.^2)./eps_silver)).*exp(-1i*x(i)*kx); % Define the Green's function
%     H(i) = sum(G(kx).*diff(kx)); % Integrate over the contour
% end
% 
% plot(x,abs(H(:,1))/MAX)