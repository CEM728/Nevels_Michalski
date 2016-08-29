clear
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
k_air = omega*sqrt(mu_0*eps_0); % propagation constant of air
k_silver = omega * sqrt(mu_0*eps_0*eps_silver); % propagation constant of silver
%%
len = 1e3; % Vector Length
%%
kxx = horzcat(linspace(0*k_air,.5*k_air,len/4),...
    linspace(.5*k_air + 1e-6, k_air,len/4),...
    linspace(k_air + 1e-6, 1.5*k_air,len/4),...
    linspace(1.5*k_air + 1e-6, 1e1*k_air,len/4)); % 1e-6 used to avoid duplicates in the array

kxy = horzcat(linspace(-1e2*k_air,-.5*k_air,len/4),...
    linspace(-.5*k_air + 1e-6, -.1*k_air,len/4),...
    linspace(-.1*k_air + 1e-6, 0,len/4),...
    linspace(0 + 1e-6, 1e1*k_air,len/4));

% Find the branch point location on the kx x-axis
dif = abs(kxx - k_air);
bp1_loc = find(dif == min(dif)); % Index in kxx with the nearest value of k_air

x = horzcat(linspace(1e-2*lambda,1e0*lambda,len),...
    linspace(1e0*lambda,1e4*lambda,len));
H_c1 = zeros(length(x),1); % Initialize the magnetic field vector
H_c2 = zeros(length(x),1); % Initialize the magnetic field vector
%% Define Contour

%               | k1 |
%               |    |
%             C1|    |c2
%               |    |
%       ------> |    V <--------
%       bottom  |    | top sheet
%               ^    |
% Define real and imaginary start points for c1
%
% Path location in terms of re_kx in kx plane
c1_start_real = kxx(bp1_loc);
c1_start_imag = -kxx(bp1_loc)*1e1;
diff = abs(kxy - c1_start_imag);
c1_sty_loc = find(diff == min(diff)); % Starting point on the Im_kx axis

%
% Define real and imaginary end points for c1
%
c1_end_real = kxx(bp1_loc);
diff = abs(kxy - 0);
c1_eny_loc = find(diff == min(diff)); % Ending point on the Im_kx axis
% c1_end_imag = kxy(c1_eny_loc) - 1e-6;
c1_end_imag = 0;
%
% Make C1 Contour
c1_real = c1_start_real*ones(len*1e2,1);
c1_imag = linspace(c1_start_imag, c1_end_imag,len*1e2);
%
c1 = horzcat(c1_real, c1_imag'); % concatenate real and imaginary parts
kx_c1 = c1_real(:,1) + 1i*c1_imag(1,:)';


% Define real and imaginary start points for c2
%
c2_start_real = kxx(bp1_loc);
c2_start_imag = 0;
%
% Define real and imaginary end points for c2
%
c2_end_real = kxx(bp1_loc);
c2_end_imag = -kxx(bp1_loc)*1e1;
%
% Make c2 Contour
c2_real = c2_start_real*ones(len*1e2,1);

c2_imag = linspace(c2_start_imag, c2_end_imag,len*1e2);
c2 = horzcat(c2_real, c2_imag');
kx_c2 = c2_real(:,1) + 1i*c2_imag(1,:)';


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
% Branch Cut Curve
hyp_silver = imag(k_silver^2)./(2*kxx); % Hyperbolic cruve for silver

% Intersection of Branch cut with vertical cut
y_int = imag(k_silver^2)/(2*k_air);

%% Display Progress Bar
h = waitbar(0,'Please wait...');
steps = length (x);

%% Integrate
for i = 1 : length (x)
    % C1 lies totally on the bottom sheet of k1
    % C1 partially lies on the bottom sheet of k2 ( until y_int)
    %    and partially on the top sheet
    
    % Integrate on left edge
   
    for j = 1 : length(kx_c1)
        
        % Enfore bottom sheet on k2
        kz1_c1(j) = sqrt(k_air^2 - kx_c1(j)^2);
        
        if abs(imag(kx_c1(j))) < abs(y_int)
            
            % Enfore bottom sheet on k2
            kz2_c1(j) = -1i*sqrt(kx_c1(j)^2 - k_silver^2);

        else
            % Top sheet of kz2
            kz2_c1(j) = sqrt(k_silver^2 - kx_c1(j)^2);

        end
    end

    D_c1 = D(kz1_c1, kz2_c1);
    G_c1 = 1./D_c1;
    integrand_1 = G_c1.*exp(-1i*kx_c1*x(i));
    H_c1(i) = trapz(integrand_1);   
    
    % Right Edge
    for j = 1 : length(kx_c2)
        
        % Enfore top sheet on k1
        kz1_c2(j) = -1i*sqrt(kx_c2(j)^2 - k_air^2);

        if abs(imag(kx_c2(j))) < abs(y_int)
            
            % Enfore top sheet on k2
            kz2_c1(j) = sqrt(k_silver^2 - kx_c1(j)^2);

        else
            
            % Bottom sheet of kz2
            kz2_c1(j) = -1i*sqrt(kx_c1(j)^2 - k_silver^2);

        end
    end
       
    D_c2 = D(kz1_c2, kz2_c2);
    G_c2 = 1./D_c2;
    integrand_2 = G_c2.*exp(-1i*kx_c2*x(i));
    H_c2(i) = trapz(integrand_2);
    
    % computations take place here
    waitbar(i / steps);
    
end

close(h)
%% Define H by EQ. 1 in [1]

H = -k_air/(2*pi*eta_0)*(H_c1 - H_c2); % C2 Contour must be zero anyway, therefore taken out

%% Plot Figure

figure('Name','Decay of the Creeping Wave',...
    'Position', [876   214   630   641]); % Size according to the paper

loglog(x/lambda,abs(H),'LineWidth',1.4,'Color','black')

% Set figure background to white
set(gcf,'Color','white');

% Create xlabel
ylabel('$\vert Creeping ~Wave\vert$',...
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

grid on

% Set y-limits according to the paper
ylim([10e-6 10e0])

title('Decay Plot of Creeping Wave part');

%% Save as a tikZ object

cleanfigure();
matlab2tikz('filename',sprintf('figures/kx_int_633.tex'),'showInfo', false)

%% Save data for use in other files

save data/kx_int_633.mat % Save data to plot the branch cuts
save('data/fields_633.mat', 'H', 'x')              % Save H and x for plotting along SPP

%% Plot the integrands

% figure(2)
% integrand1 = G_c1.*exp(-1i*1e2*kx_c1);
% loglog(imag(kx_c1), real(integrand1),'LineWidth',1.4,'Color','black')
% set(gcf,'Color','white');
% % hold on
% % loglog(imag(kx_c1), imag(integrand1),'LineWidth',1.4,'Color','black',...
% %     'linestyle', '--');
% ylabel('Integrand along C1',...
%     'HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'FontSize',12,...
%     'Interpreter','latex');
% 
% % Create xlabel
% xlabel('$\Im k_x$',...
%     'HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'FontSize',12,...
%     'Interpreter','latex');
% grid on
% % legend('Real Part', 'Imaginary Part')
% 
% 
% figure(2)
% integrand2 = G_c2.*exp(-1i*1e0*kx_c2);
% loglog(imag(kx_c2), real(integrand2),'LineWidth',1.4,'Color','black')
% % hold on
% % loglog(imag(kx_c2), imag(integrand2),'LineWidth',1.4,'Color','black',...
% %     'linestyle', '--');
% set(gcf,'Color','white');
% 
% ylabel('Integrand along C2',...
%     'HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'FontSize',12,...
%     'Interpreter','latex');
% 
% % Create xlabel
% xlabel('$\Im k_x$',...
%     'HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'FontSize',12,...
%     'Interpreter','latex');
% % legend('Real Part', 'Imaginary Part')
% grid on

%%
toc