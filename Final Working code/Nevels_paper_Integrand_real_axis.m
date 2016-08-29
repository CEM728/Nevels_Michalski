clear all
close all;tic

lambda = 1;%633e-9; % Red light wavelength
eps_silver =  -18.295 - 1i*0.48085; % Johnson & Christy,1972 (refractiveindex.info) at 633 nm
% eps_silver = -265.06 - 1i*29.436;  % @ 2500 nm
load em_constants.mat % Contains varepsilon, mu and c
eps_0 = epsilon_0;
eps_1 = 1 ;
c = 1/sqrt(mu_0*eps_0);
omega = 2*pi*c/lambda; % angular frequency
eta_0 = sqrt(mu_0/eps_0);
k_air = omega*sqrt(mu_0*eps_0)+ .0001*1i; % propagation constant of air
k_silver = omega * sqrt(mu_0*eps_0*eps_silver); % propagation constant of silver
kxp = k_air*sqrt(1*eps_silver/(eps_silver+ 1)); % SPP pole location
%%
len = 5e3; % Vector Length

% kx_c1 = (-20:1e-3:20)';
% kx_c2 = (-20:1e-3:20)';
kx_c1 = (1:1e-3:10)'; % Zoomed in
kx_c2 = (1:1e-3:10)';


%% Define Green's function


kz_1 = @(kx) sqrt(k_air^2 - kx.^2);
kz_2 = @(kx) sqrt(k_silver^2 - kx.^2);
D = @(kz_1, kz_2) kz_2/eps_silver + kz_1/eps_1;
% G = @(kz_1, kz_2) eps_0./D;

% On Contour C1
kz1_c1 = kz_1(kx_c1);
kz2_c1 = kz_2(kx_c1);


% On Contour c2
kz1_c2 = kz_1(kx_c2);
kz2_c2 = kz_2(kx_c2);

%%

D_c1 = D(kz1_c1, kz2_c1);
G_c1 = 1./D_c1;

D_c2 = D(kz1_c2, kz2_c2);
G_c2 = 1./D_c2;

%% Plot the integrands

figure(1)
integrand1 = G_c1;
plot(real(kx_c1), real(integrand1),'LineWidth',1.4)
hold on
plot(real(kx_c1), imag(integrand1),'LineWidth',1.4)
set(gcf,'Color','white');
% Plot branch points
plot(real(k_air),imag(k_air),'Marker','o',...
                'LineWidth',1.6,...
                'MarkerEdgeColor','black',...
                'MarkerSize',5)
plot(real(kxp),0,'Marker','x',...
                'LineWidth',1.6,...
                'MarkerEdgeColor','green',...
                'MarkerSize',5)

ylabel('Integrand along Real Axis',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create xlabel
xlabel('$\mathbf{k_x}$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create Legend
legend({'Real Part', 'Imaginary Part',...
    'Branch Point','Pole'},...
    'FontSize',10,...
    'Interpreter','latex');

% Create Annotation box
str = '$\frac{x}{\lambda_0} = 1\times 10^1$';
dim = [.2 .2 .2 .2];
annotation('textbox',dim,'String',str,...
    'FitBoxToText','on',...
    'FontSize',10,...
    'Interpreter','latex');
grid on
% cleanfigure();
% matlab2tikz('filename',sprintf('figures/real_axis_integrand.tex'),'showInfo', false)
%%
toc