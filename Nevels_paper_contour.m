% load nevels_michalski_creeping_wave_rev7.mat
clear all
load data/kx_int_633.mat
close all; 
% SPP pole location
kxp = k_air*sqrt(1*eps_silver/(eps_silver+ 1));

% On Contour C1
kz1_c1 = kz_1(kx_c1);
kz2_c1 = kz_2(kx_c1);


% On Contour c2
kz1_c2 = kz_1(kx_c2);
kz2_c2 = kz_2(kx_c2);

% Branch Cut Curves
yy_1 = imag(k_silver^2)./(2*kxx); % Silver
yy_2 = imag(k_air^2)./(2*kxx); % Air


%% Plot pole, branch points and contours

% Plot pole
figure('Name','Contour');
hold on 
plot(real(kxp),imag(kxp),'Marker','x',...
                'Color','red',...
                'LineWidth',2,...
                'MarkerEdgeColor','red',...
                'MarkerSize',5)

set(gcf,'Color','white');
box on
% Plot Contour
plot(real(kx_c1),imag(kx_c1),'LineWidth',1.1,'Color','black') % Left edge 
% plot(real(kx_c2),imag(kx_c2),'LineWidth',1.1,'Color','black') % Right edge

% Plot branch points
plot(real(k_air),imag(k_air),'Marker','o',...
                'LineWidth',1.6,...
                'MarkerEdgeColor','black',...
                'MarkerSize',5)
plot(real(k_silver),imag(k_silver),'Marker','o',...
                'LineWidth',1.6,...
                'MarkerEdgeColor','green',...
                'MarkerSize',5)
            
ylabel('$\Im(k_x)$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');

% Create xlabel
xlabel('$\Re(k_x)$',...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',12,...
    'Interpreter','latex');  


str = '$\lambda_0 = 633~nm$';
dim = [.2 .2 .2 .2];
annotation('textbox',dim,'String',str,...
    'FitBoxToText','on',...
    'FontSize',10,...
    'Interpreter','latex');
% Plot Sommerfeld Branch Cut
plot(kxx, yy_1,'LineWidth',1.4,'linestyle', '-.')
plot(kxx, yy_2,'LineWidth',1.4,'linestyle', ':')

legend({'Pole', 'Vertical Cut', '$BP_1$','$BP_2$', 'Sommerfeld $BC_1$',...
        'Sommerfeld $BC_2$'},...
    'FontSize',10,...
    'Interpreter','latex',...
    'Location','SouthEast');
xlim([ 0 7])
ylim([ -30 .5])
% Save Figure
cleanfigure();
matlab2tikz('filename',sprintf('figures/contour_633.tex'),'showInfo', false)