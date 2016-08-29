%% This program performs the integration in kx plane with function handles
clear
close all;tic
%% EM Paramters
lambda = 1;%633e-9; % Red light wavelength
eps_silver =  -18.295 - 1i*0.48085; % Johnson & Christy,1972 (refractiveindex.info) at 633 nm
load em_constants.mat % Contains varepsilon, mu and c
eps_0 = epsilon_0;
eps_1 = 1 ;
c = 1/sqrt(mu_0*eps_0);
omega = 2*pi*c/lambda; % angular frequency
eta_0 = sqrt(mu_0/eps_0);
k_air = omega*sqrt(mu_0*eps_0); % propagation constant of air
k_silver = omega * sqrt(mu_0*eps_0*eps_silver); % propagation constant of silver
%%
len = 5e1;
x = horzcat(linspace(1e-2*lambda,1e0*lambda,len/2),...
    linspace(1e0*lambda,1e4*lambda,len/2));

kz_1 = @(kx) sqrt(k_air^2 - kx.^2);
kz_2 = @(kx) sqrt(k_silver^2 - kx.^2);
D = @(kx) kz_2(kx)/eps_silver + kz_1(kx)/eps_1;
G = @(kx) 1./D(kx);
Integrand  = @(kx,x) G(kx).*exp(1i*kx*x);
for i = 1 : len
    H1(i) = integral(@(kx) Integrand ( kx, x(i)), k_air-100*1i,k_air-50*1i,...
        'Waypoints',[k_air,0],...
        'RelTol',0,'AbsTol',1e-12);
%     H2(i) = integral(@(kx) Integrand ( kx, x(i)), k_air-1e-6*1i,k_air-100*k_air*1i,...
%         'Waypoints',[k_air,0],...
%         'RelTol',0,'AbsTol',1e-12);
end
loglog(x/lambda,abs(H1))