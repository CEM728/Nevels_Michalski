% clear all
close all

% lambda = 633e-9; % Green light wavelength
lambda = 2500e-9;
% eps_silver =  -18.295 - 1i*0.48085; % Johnson & Christy,1972 (refractiveindex.info) at 633 nm
eps_silver = -265.06 - 1i*29.436;  % @ 2500 nm 
load em_constants.mat % Contains varepsilon, mu and c
eps_0 = epsilon_0;
omega = 2*pi*c/lambda; % angular frequency

len = 1e4;
k_air = 2*pi/lambda; % propagation constant of air
k_silver = omega * sqrt(mu_0*epsilon_0*eps_silver); % propagation constant of silver
%
kxx = linspace(0*k_air,1e4*k_air,1e2*len);
kxy = linspace(-1e3*k_air,1e1*k_air,1e2*len);    
% kx = meshgrid(kxx,kxy); % [real(kx),imag(kx)]
dif = abs(kxx - k_air);
pole_location = find(dif == min(dif)); % Index in kxx with the nearest value of k_air

x1 = linspace(1e-2*lambda,1e4*lambda,len);

load test_rev5.mat
% load test_2500nm.mat
% x = linspace(1e-2*lambda,1e4*lambda,1e1*len); % Space vector definition
% H = zeros(length(x)); % Initialize the magnetic field vector

%%
kxp = k_air*sqrt(1*eps_silver/(eps_silver+ 1)); % SPP pole location
% beta defintions
kz_1 = @(kx) sqrt(k_air^2 - kx.^2);
kz_2 = @(kx) sqrt(k_silver^2 - kx.^2);
%
D_tilde = @(kx) -kx.*(1./(eps_silver*kz_2(kx)) + 1./(1*kz_1(kx)));
R_p = 1./D_tilde(kxp); 
SPP = -2*pi*1i*R_p*exp(-1i*kxp*x1);

%%
loglog(x1/lambda,abs(SPP),'LineWidth',1.4,'Color','black',...
    'LineStyle' , '--');
hold on
loglog(x,abs(H),'LineWidth',1.4,'Color','black');
set(gcf,'Color','white');
hold off
ylabel('$\vert Plasmon Wave\vert$',...    
'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create xlabel
xlabel('$\frac{x}{\lambda}$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');  

title('Decay Plot of Plasmon wave Wave part');
legend('Plasmon Wave', 'Creeping Wave');
ylim([10e-6 10e1])
xlim([1e-2 1e4])
%%
cleanfigure();
matlab2tikz('filename',sprintf('nevels_michalski_both_plots.tex'))