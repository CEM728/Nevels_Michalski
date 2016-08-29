%% Transform the integral in s-plane

clear
close all;clf;clc;tic
%%
lambda = 1;% 9000e-9; % Red light wavelength
% eps_silver =  -18.295 - 1i*0.48085; % Johnson & Christy,1972 (refractiveindex.info) at 633 nm
eps_silver = -265.06 - 1i*29.436;  % @ 2500 nm
% eps_silver = -3110.8 - 1i*1116.5; %@9000 nm
load em_constants.mat % Contains varepsilon, mu and c
eps_0 = epsilon_0;
eps_1 = 1;
c = 1/sqrt(mu_0*eps_0);
omega = 2*pi*c/lambda; % angular frequency
eta_0 = sqrt(mu_0/eps_0);
k_air = 2*pi/lambda; % propagation constant of air
k_silver = omega * sqrt(mu_0*epsilon_0*eps_silver); % propagation constant of silver
kxp = k_air*sqrt(1*eps_silver/(eps_silver+ 1)); % SPP pole location
sxp = exp(1i*pi/4)*sqrt(kxp - k_air);

kx = @(s) k_air - 1i*s.^2;

%% Array Sizes
len = 1e4; % Vector Length
% %%
% x = horzcat(linspace(1e-2*lambda,1e0*lambda,len/2),...
%     linspace(1e0*lambda,1e4*lambda,len/2));
x = 5*lambda;
% Initialize
integrand = zeros(1,len);
Creep = zeros(1,len);

% beta defintions
kz_1 = @(kx) -sqrt(k_air^2 - kx.^2);
kz_2 = @(kx) sqrt(k_silver^2 - kx.^2);
%
D_tilde = @(kx) -kx.*(1./(eps_silver*kz_2(kx)) + 1./(1*kz_1(kx)));
R_p = 1./D_tilde(kxp); 

kz1 = @(s) -s.*sqrt(1i*2*k_air + s.^2);
kz2 = @(s) -sqrt(k_silver^2 - k_air^2 + 1i*2*s.^2*k_air + s.^4);
D_tilde1 = @(s) -(k_air - 1i*s.^2).*(1./(eps_silver*kz2(s)) + 1./(1*kz1(s)));
R_p1 = 1./D_tilde1(sxp);
B_p = 1j*R_p1./sxp;

% kz1 = kz1(s);
% kz2 = kz2(s);

% Integrands

% With pole
F = @(s) -2.*kz1(s)/eps_1./((kz2(s)/eps_silver).^2 - (kz1(s)/eps_1).^2);

% Without pole
F_p = @(s) (F(s)- B_p.*s./(s.^2 - sxp^2));

%
I1 = @(s) F(s).*exp(-x.*s.^2).*s;
I2 = @(s) (F(s)- B_p.*s./(s.^2 - sxp^2)).*exp(-x.*s.^2).*s;
Ip = @(s) B_p./(s.^2 - sxp^2).*exp(-x.*s.^2).*s.^2;
S = linspace(0,2,len);
F1 = I1(S);
F2 = I2(S);
% F2 = Ip(S);

close(1)
%%
figure('Name','Pole Effect',...
    'Position', [876   214   630   641]); % Size according to the paper
hold on
plot(S,real(F1),'LineWidth',1.1,'Color','black')
plot(S,imag(F1),'LineWidth',1.1,'Color','black','linestyle', '--')
% loglog(x/lambda,abs(Creep),'LineWidth',1.4,'Color','black')
% loglog(x, abs(H)/abs(max(H)),'LineWidth',1.4,'Color','black')
set(gcf,'Color','white');

ylabel('Integrand Amplitude',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create xlabel
xlabel('$s$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% ylim([10e-6 10e1])
% xlim([1e-2 1e4])
title('Creeping Wave Integrand with the pole');

cleanfigure();
matlab2tikz('filename',sprintf('nevels_michalski_decay_plot_S_trans_pole_effect_9000.tex'),'showInfo', false)


%%
figure('Name','Pole Effect Substracted',...
    'Position', [876   214   630   641]); % Size according to the paper
hold on
plot(S,real(F2),'LineWidth',1.1,'Color','black')
plot(S,imag(F2),'LineWidth',1.1,'Color','black','linestyle', '--')
% loglog(x/lambda,abs(Creep),'LineWidth',1.4,'Color','black')
% loglog(x, abs(H)/abs(max(H)),'LineWidth',1.4,'Color','black')
set(gcf,'Color','white');

ylabel('Integrand Amplitude',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create xlabel
xlabel('$s$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% ylim([10e-6 10e1])
% xlim([1e-2 1e4])
title('Creeping Wave Integrand with the pole');
cleanfigure();
matlab2tikz('filename',sprintf('nevels_michalski_decay_plot_pole_effect_9000_subtracted_.tex'),'showInfo', false)
% %%
% 
% save nevels_michalski_creeping_S_transform.mat % Save data to plot the branch cuts
save('test_S_transform_2500.mat', 'Creep', 'x')              % Save H and x for plotting along SPP
toc