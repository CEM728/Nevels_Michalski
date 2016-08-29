%% Transform the integral in s-plane

clear
close all;clf;clc;tic
%%
lambda = 2500e-9; % Red light wavelength
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


%% Array Sizes
len = 5e3; % Vector Length
% %%
x = horzcat(linspace(1e-2*lambda,1e0*lambda,len/2),...
    linspace(1e0*lambda,1e4*lambda,len/2));
% x = 5*lambda;
% S = linspace(0, sxp*1e2, len);
% Initialize
integrand = zeros(1,len);
Creep = zeros(1,len);

kz1 = @(s) -s.*sqrt(1i*2*k_air + s.^2);
kz2 = @(s) sqrt(k_silver^2 - k_air^2 + 1i*2*s.^2*k_air + s.^4);

% kz1 = kz1(s);
% kz2 = kz2(s);

F = @(s) -2.*kz1(s)/eps_1./((kz2(s)/eps_silver).^2 - (kz1(s)/eps_1).^2);
I1 = @(x,s) F(s)*exp(-x*s.^2).*s;
S = linspace(0,10,len);
F1 = F(S);

S1 = linspace(0,15,len);
h = waitbar(0,'Processing...');
steps = len;
for i = 1 : len
    for j = 1 : len
        I1 = @(s) F(s)*exp(-x(i)*s.^2).*s;
        integrand(j) = I1(S1(j));
    end
    Creep(i) = -2*1i*exp(-1i*x(i)*k_air).*trapz(integrand);
    % computations take place here
    waitbar(i / steps);
    % end
    
end
close(h)
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
matlab2tikz('filename',sprintf('nevels_michalski_decay_plot_S_trans_pole_effect.tex'),'showInfo', false)

%%
figure('Name','Creeping Wave Decay',...
    'Position', [876   214   630   641]); % Size according to the paper
loglog(x/lambda,abs(Creep),'LineWidth',1.4,'Color','black')
set(gcf,'Color','white');

ylabel('Creeping Wave',...
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

% ylim([10e-6 10e1])
xlim([1e-2 1e4])
title('Creeping Wave Decay');
cleanfigure();
matlab2tikz('filename',sprintf('nevels_michalski_decay_plot_S_trans_2500.tex'),'showInfo', false)
%%

save nevels_michalski_creeping_S_transform.mat % Save data to plot the branch cuts
save('test_S_transform_2500.mat', 'Creep', 'x')              % Save H and x for plotting along SPP
toc