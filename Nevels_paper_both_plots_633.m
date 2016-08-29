% clear all
close all
%%
LW = 'linewidth'; lw = 1.6;

lambda = 633e-9; % Green light wavelength
eps_silver =  -18.295 - 1i*0.48085; % Johnson & Christy,1972 (refractiveindex.info) at 633 nm
% eps_silver = -265.06 - 1i*29.436;  % @ 2500 nm
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

load data/fields.mat
% load data/fields_2500.mat
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
loglog(x1/lambda,abs(SPP),LW,lw,...
     'LineStyle' , '-.');
hold on
loglog(x,2e-1*abs(H),LW,lw);
set(gcf,'Color','white');
ax = gca;
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
hold off
ylabel('\textbf{\textit{Amplitude}}',...    
'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create xlabel
xlabel('$\mathbf{\mathit{\frac{x}{\lambda}}}$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');  

% title('Logarithmic Decay Plot Comparison of Creeping and Plasmon wave',...
%     'FontWeight','bold',...
%     'FontSize',12,...
%     'Interpreter','latex');
legend({'Plasmon Wave', 'Creeping Wave'},...
    'FontSize',10,...
    'Interpreter','latex');
str = '$\lambda_0 = 633~nm$'
dim = [.2 .2 .2 .2];
annotation('textbox',dim,'String',str,...
    'FitBoxToText','on',...
    'FontSize',10,...
    'Interpreter','latex');

ylim([10e-6 10e0])
xlim([1e-2 1e4])
%%
cleanfigure();
matlab2tikz('filename',sprintf('figures/kx_int_overlay_633.tex'))