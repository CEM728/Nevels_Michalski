load nevels_michalski_creeping_wave_2500nm.mat
close all;clf
% SPP pole location
kxp = k_air*sqrt(1*eps_silver/(eps_silver+ 1));

% On Contour C1
kz1_c1 = kz_1(kx_c1);
kz2_c1 = kz_2(kx_c1);

% On Contour C2
kz1_c2 = kz_1(kx_c2);
kz2_c2 = kz_2(kx_c2);

% On Contour C3
kz1_c3 = kz_1(kx_c3);
kz2_c3 = kz_2(kx_c3);

% Branch Cut Curves
yy_1 = imag(k_silver^2)./(2*kxx); % 


%% Plot pole, branch points and contours

% Plot pole
figure('Name','Contour');
plot(real(kxp),imag(kxp),'mo',...
                'LineWidth',2,...
                'MarkerEdgeColor','none',...
                'MarkerFaceColor',[0.5 0.5 0.5],...
                'MarkerSize',5)
hold on 

% Plot Contour
plot(real(kx_c1),imag(kx_c1),'LineWidth',1.1,'Color','black') % Left edge 
plot(real(kx_c3),imag(kx_c3),'LineWidth',1.1,'Color','black') % Right edge
plot(real(kx_c2),imag(kx_c2),'LineWidth',1.1,'Color','black') % Top edge

% Plot branch points
plot(real(k_air),imag(k_air),'go',...
                'LineWidth',2,...
                'MarkerEdgeColor','none',...
                'MarkerFaceColor',[0.5 0.5 0.5],...
                'MarkerSize',5)
plot(real(k_silver),imag(k_silver),'ko',...
                'LineWidth',2,...
                'MarkerEdgeColor','none',...
                'MarkerFaceColor',[0.5 0.5 0.5],...
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

% Plot Sommerfeld Branch Cut
plot(kxx, yy_1,'LineWidth',1.1,'Color','black','linestyle', '--')
% xlim([ -1 10])
% ylim([ -40 5])
% Save Figure
matlab2tikz('filename',sprintf('nevels_michalski_contour.tex'))