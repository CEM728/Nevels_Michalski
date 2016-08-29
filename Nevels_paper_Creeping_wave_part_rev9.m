clear
close all;clc
%%
lambda = 633e-9; % Red light wavelength
eps_silver =  -18.295 - 1i*0.48085; % Johnson & Christy,1972 (refractiveindex.info) at 633 nm
% eps_silver = -265.06 - 1i*29.436;  % @ 2500 nm
load em_constants.mat % Contains varepsilon, mu and c
eps_0 = epsilon_0;
eps_1 = 1;
c = 1/sqrt(mu_0*eps_0);
omega = 2*pi*c/lambda; % angular frequency
eta_0 = sqrt(mu_0/eps_0);
k_air = 2*pi/lambda; % propagation constant of air
k_silver = omega * sqrt(mu_0*epsilon_0*eps_silver); % propagation constant of silver
%% Array Sizes
len = 2e2; % Vector Length
%%
kxx = horzcat(linspace(0*k_air,.5*k_air,len/4),...
    linspace(.5*k_air + 1e-6, k_air,len/4),...
    linspace(k_air + 1e-6, 1.5*k_air,len/4),...
    linspace(1.5*k_air + 1e-6, 1e3*k_air,len/4));

kxy = horzcat(linspace(-1e2*k_air,-1e1*k_air,len/4),...
    linspace(-1e1*k_air + 1e-6, -1e0*k_air,len/4),...
    linspace(-k_air + 1e-6, -1e-1*k_air,len/4),...
    linspace(-1e-1*k_air + 1e-6, 1e2*k_air,len/4));


% kxx = linspace(-1e1*k_air,1e1*k_air,len);
% kxy = linspace(-1e4*k_air,1e3*k_air,len);

x = horzcat(linspace(1e-2*lambda,1e0*lambda,len/2),...
    linspace(1e0*lambda,1e4*lambda,len/2));

% x = linspace(1e-2*lambda,1e4*lambda,len);

H_c1 = zeros(len,1); % Initialize the magnetic field vector
H_c2 = zeros(len,1); % Initialize the magnetic field vector

%% Find the branch point location on the kx x-axis
dif = abs(kxx - k_air);
bp1_loc = find(dif == min(dif)); % Index in kxx with the nearest value of k_air

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
c1_start_imag = -kxx(bp1_loc)*1e3;
diff = abs(kxy - c1_start_imag);
c1_sty_loc = find(diff == min(diff)); % Starting point on the Im_kx axis

%
% Define real and imaginary end points for c1
%
c1_end_real = kxx(bp1_loc);
diff = abs(kxy - 0);
c1_eny_loc = find(diff == min(diff)); % Ending point on the Im_kx axis
% c1_end_imag = kxy(c1_eny_loc)/len^2;
c1_end_imag = 1e-6;
%
% Make C1 Contour
c1_real = c1_start_real*ones(len,1);
c1_imag = linspace(c1_start_imag, c1_end_imag,len);
%
c1 = horzcat(c1_real, c1_imag'); % concatenate real and imaginary parts
kx_c1 = c1_real(:,1) + 1i*c1_imag(1,:)';


% Define real and imaginary start points for c2
%
c2_start_real = kxx(bp1_loc);
% c2_start_imag = kxy(c1_eny_loc)/len^2;
c2_start_imag = 1e-6;
%
% Define real and imaginary end points for c2
%
c2_end_real = kxx(bp1_loc);
c2_end_imag = -kxx(bp1_loc)*1e3;
%
% Make c2 Contour
c2_real = c2_start_real*ones(len,1);
c2_imag = linspace(c2_start_imag, c2_end_imag,len);
c2 = horzcat(c2_real, c2_imag');
kx_c2 = c2_real(:,1) + 1i*c2_imag(1,:)';


%% Define Green's function


kz_1 = @(kx) sqrt(k_air^2 - kx.^2);
kz_2 = @(kx) sqrt(k_silver^2 - kx.^2);
D = @(kz_1, kz_2) kz_2/eps_silver + kz_1/eps_1;
G = @(kz_1, kz_2) 1./D;

% On Contour C1
kz1_c1 = kz_1(kx_c1);
kz2_c1 = kz_2(kx_c1);


% On Contour c2
kz1_c2 = kz_1(kx_c2);
kz2_c2 = kz_2(kx_c2);

%%
% Branch Cut Curve
hyp_silver = imag(k_silver^2)./(2*kxx); % Hyperbolic curve for silver

% Intersection on C1
[X0,Y0] = intersections(real(kx_c1),imag(kx_c1),kxx,hyp_silver,1);
Silver_branch_cut_loc_c1 = Y0;

% Intersection on c2
[X1,Y1] = intersections(real(kx_c2),imag(kx_c2),kxx,hyp_silver,1);
Silver_branch_cut_loc_c2 = Y1;

%% Display Progress Bar
h = waitbar(0,'Processing...');
steps = length (x);

%% Enforce proper square-root terms
for i = 1 : len
    % C1 lies totally on the bottom sheet of k1
    % C1 partially lies on the bottom sheet of k2 ( until Silver_bc_l_c1)
    %    and partially on the top sheet
    
    % For k1
    %   Re(kz1) > 0, Im(kz1) > 0 ---- Bottom Sheet
    % For k2
    %   Re(kz2) > 0, Im(kz2) > 0 for kx < Silver_loc_c1 ---- Bottom Sheet
    %   Re(kz2) > 0, Im(kz2) < 0 for kx > Silver_loc_c1 ---- Top Sheet
    
    % Left Path
    for j = 1 : length(kx_c1)
        
        % Bottom Sheet on kz1
        % Im(kz1) > 0, Re(kz1) > 0
        if real(kz1_c1(j)) < 0
            kz1_c1(j) = -conj(kz1_c1(j));
        end
        if imag(kz1_c1(j)) < 0
            kz1_c1(j) = conj(kz1_c1(j));
        end
        
        if abs(imag(kx_c1(j))) < abs(Silver_branch_cut_loc_c1)
            % Enfore bottom sheet on k2
            % Im(kz2) > 0, Re(kz2) > 0
            
            if real(kz2_c1(j)) < 0
                kz2_c1(j) = -conj(kz2_c1(j));
            end
            
            if imag(kz2_c1(j)) < 0
                kz2_c1(j) = conj(kz2_c1(j));
            end
        else
            % Top sheet of kz2
            % Im(kz2) < 0, Re(kz2) > 0
            
            if real(kz2_c1(j)) < 0
                kz2_c1(j) = -conj(kz2_c1(j));
            end
            
            if imag(kz2_c1(j)) > 0
                kz2_c1(j) = conj(kz2_c1(j));
            end
        end
        
        
    end
    
    % Right Path
    % c2 lies totally on the top sheet of k1
    % c2 partially lies on the top sheet of k2 ( until Silver_bc_l_c2)
    %    and partially on the bottom sheet
    
    % For k1
    %   Re(kz1) < 0, Im(kz1) < 0 ---- Top Sheet
    % For k2
    %   Re(kz2) > 0, Im(kz2) < 0 for kx < Silver_loc_c2 ---- Top Sheet
    %   Re(kz2) > 0, Im(kz2) > 0 for kx > Silver_loc_c2 ---- Bottom Sheet
    for j = 1 : length(kx_c2)
        
        % Enfore top sheet on k1
        % Im(kz1) < 0, Re(kz1) < 0
        if real(kz1_c2(j)) > 0
            kz1_c2(j) = -conj(kz1_c2(j));
        end
        
        if imag(kz1_c2(j)) > 0
            kz1_c2(j) = conj(kz1_c2(j));
        end
        
        
        if abs(imag(kx_c2(j))) < abs(Silver_branch_cut_loc_c2)
            % Enfore top sheet on k2
            % Im(kz2) < 0, Re(kz2) > 0
            if real(kz2_c2(j)) < 0
                kz2_c2(j) = -conj(kz2_c2(j));
            end
            
            if imag(kz2_c2(j)) > 0
                kz2_c2(j) = conj(kz2_c2(j));
            end
        else
            % Bottom sheet of kz2
            % Im(kz2) > 0, Re(kz2) > 0
            
            if real(kz2_c2(j)) > 0
                kz2_c2(j) = -conj(kz2_c2(j));
            end
            
            if imag(kz2_c2(j)) < 0
                kz2_c2(j) = conj(kz2_c2(j));
            end
        end
    end
    
    
    % Left Path
    D_c1 = D(kz1_c1, kz2_c1);
    G_c1 = 1./D_c1;
    
    % Right Path
    D_c2 = D(kz1_c2, kz2_c2);
    G_c2 = 1./D_c2;
    
    
    
    integrand_1 = G_c1.*exp(-1i*k_air*kx_c1*x(i));
    H_c1(i) = trapz(integrand_1);
    integrand_3 = G_c2.*exp(-1i*k_air*kx_c2*x(i));
    H_c2(i) = trapz(integrand_3);
end
% computations take place here
waitbar(i / steps);
% end
close(h)
%% Define H by EQ. 1 in [1]
H = -k_air/(2*pi*eta_0)*(H_c1  +  H_c2); % C2 Contour must be zero anyway, therefore taken out
%% Plot Figure
figure('Name','Decay of the Creeping Wave',...
    'Position', [876   214   630   641]); % Size according to the paper

loglog(x/lambda, abs(H)/abs(max(H)),'LineWidth',1.4,'Color','black')
% loglog(x, abs(H)/abs(max(H)),'LineWidth',1.4,'Color','black')
set(gcf,'Color','white');

ylabel('$\vert Creeping Wave\vert$',...
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

ylim([10e-5 10e1])

title('Decay Plot of Creeping Wave part');

cleanfigure();
matlab2tikz('filename',sprintf('nevels_michalski_decay_plot_rev7.tex'),'showInfo', false)
%%
save nevels_michalski_creeping_wave_rev7.mat % Save data to plot the branch cuts
save('test_rev5.mat', 'H', 'x')              % Save H and x for plotting along SPP
%% Transpose All variables for all variables
% kz_1 = kz_1.';
% kz_2 = kz_2.';
% D = D.';
% G = G.';
% H = H.';
% kx = kx.';