%% Plot Riemann Surface of the Integrand of the Sommerfeld Integral
clear all; close all
% set(gcf,'Color','white'); % Set background color to white
% title('Riemann Surface of $f(z) = \sqrt{z}$','Interpreter','latex')
% set (gca,'FontName','times new roman') % Set axes fonts to Times New Roman
% ax = gca;
% xlabel('$\Re(z)$','Interpreter','latex');
% ylabel('$\Im(z)$','Interpreter','latex');
% material dull; % Set reflectivity of the surface to dull
%%
lambda = 1; %633e-9; % Red light wavelength
eps_silver =  -18.295 - 1i*0.48085; % Johnson & Christy,1972 (refractiveindex.info) at 633 nm
load em_constants.mat % Contains varepsilon, mu and c
clear c
eps_0 = epsilon_0;
c1 = 1/sqrt(mu_0*eps_0);
omega = 2*pi*c1/lambda; % angular frequency
eta_0 = sqrt(mu_0/eps_0);

k_air = 2*pi/lambda; % propagation constant of air
k_silver = omega * sqrt(mu_0*epsilon_0*eps_silver);

len = 1e2; % Vector Length
%%
% kxx = linspace(0*k_air,1e1*k_air,len);
% kxy = linspace(-1e1*k_air,1e1*k_air,len);

kxx = linspace(-10,10,len);
kxy = linspace(-10,10,len);
[kx, ky] = meshgrid(kxx,kxy);
K = kx + 1i*ky;
ex = 1;
kz_1 = @(kx) sqrt(k_air^2 - kx.^2);
kz_2 = @(kx) sqrt(k_silver^2 - kx.^2);
D = @(kz_1, kz_2) kz_2/eps_silver + kz_1/eps_0;
D1 = @(k) sinh(k);
D1 = sinh(K);
% On Contour C1
kz1_c1 = kz_1(K);
kz2_c1 = kz_2(K);
D_c1 = D(kz1_c1, kz2_c1);
G = 1./D_c1.*exp(-1i*K*ex);

[h,c] = dcolor(real(K), imag(K), D_c1, 'grid');
% shading interp
caxis( [ min(c(:)), max(c(:)) ] )
axis square
colormap viridis
figure(2)
PhasePlot(K,D_c1);
colormap viridis
% contour(real(K), imag(K), imag(D_c1))