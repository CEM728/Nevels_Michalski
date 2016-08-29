% clear
close all;clc

lambda = 1; %633e-9; % Red light wavelength
eps_silver =  -18.295 - 1i*0.48085; % Johnson & Christy,1972 (refractiveindex.info) at 633 nm
load em_constants.mat % Contains varepsilon, mu and c
eps_0 = epsilon_0;
c = 1/sqrt(mu_0*eps_0);
omega = 2*pi*c/lambda; % angular frequency
eta_0 = sqrt(mu_0/eps_0);

k_air = 2*pi/lambda; % propagation constant of air
k_silver = omega * sqrt(mu_0*epsilon_0*eps_silver); % propagation constant of silver
%%
len = 1e2; % Vector Length
%%
kxx = linspace(-1e1*k_air,1e4*k_air,1e4*len);
kxy = linspace(-1e4*k_air,1e1*k_air,1e4*len);    

% Find the branch point location on the kx x-axis    
dif = abs(kxx - k_air);
bp1_loc = find(dif == min(dif)); % Index in kxx with the nearest value of k_air

x = linspace(1e-2*lambda,1e4*lambda,1e2*len);
H_c1 = zeros(length(x),1); % Initialize the magnetic field vector
H_c2 = zeros(length(x),1); % Initialize the magnetic field vector
H_c3 = zeros(length(x),1); % Initialize the magnetic field vector
%% Define Contour
%                 C2
%               --->--
%               | k1 |
%               |    |
%             C1|    |C3
%               |    |
%       ------> |    V <--------
%       bottom  |    | top sheet
%               ^    |
% Define real and imaginary start points for c1
%
c1_loc = bp1_loc - 1; % Path location in terms of re_kx in kx plane
c1_start_real = kxx(bp1_loc - 1);
c1_start_imag = -kxx(bp1_loc - 1)*1e3;
diff = abs(kxy - c1_start_imag);
c1_sty_loc = find(diff == min(diff)); % Starting point on the Im_kx axis

%
% Define real and imaginary end points for c1
%
c1_end_real = kxx(bp1_loc - 1);
diff = abs(kxy - 0);
c1_eny_loc = find(diff == min(diff)); % Ending point on the Im_kx axis
c1_end_imag = kxy(c1_eny_loc - 1)/len;
%
% Make C1 Contour
c1_real = c1_start_real*ones(1e4*len,1);
c1_imag = linspace(c1_start_imag, c1_end_imag,1e4*len);
%
c1 = horzcat(c1_real, c1_imag'); % concatenate real and imaginary parts
kx_c1 = c1_real(:,1) + 1i*c1_imag(1,:)';

%
% Define real and imaginary start points for c2
%
c2_start_real = kxx(bp1_loc - 1);
c2_start_imag = kxy(c1_eny_loc + 1)/len;
diff = abs(kxy - c2_start_imag);
c2_y_loc = find(diff == min(diff)); % Starting point on the Im_kx axis
%
% Define real and imaginary end points for c2
%
c2_end_real = kxx(bp1_loc + 1);
c2_end_imag = kxy(c1_eny_loc + 1)/len;
%
%
% Make C3 Contour
%
c2_real = linspace(c2_start_real, c2_end_real,3);
c2_imag = c2_end_imag*ones(3,1);

c2 = horzcat(c2_real', c2_imag);
kx_c2 = c2_real(1,:)' + 1i*c2_imag(1,:);
% Define real and imaginary start points for c3
%
c3_start_real = kxx(bp1_loc + 1);
c3_start_imag = kxy(c1_eny_loc - 1)/len;
%
% Define real and imaginary end points for c3
%
c3_end_real = kxx(bp1_loc + 1);
c3_end_imag = -kxx(bp1_loc - 1)*1e3;
%
% Make C3 Contour
c3_real = c3_start_real*ones(1e4*len,1);
c3_imag = linspace(c3_start_imag, c3_end_imag,1e4*len);

c3 = horzcat(c3_real, c3_imag');
kx_c3 = c3_real(:,1) + 1i*c3_imag(1,:)';

%%
% Branch Cut Curve
hyp_silver = imag(k_silver^2)./(2*kxx); % Hyperbolic cruve for silver

% Intersection on C1
flip_c1_imag = flip(c1_imag);
diff1 = abs(hyp_silver - flip_c1_imag);
c1_hyp_silver_int_loc = find(diff1 == min(diff1))
Silver_branch_cut_loc_c1 = flip_c1_imag(c1_hyp_silver_int_loc)

% Intersection on C3
diff2 = abs(hyp_silver - c3_imag);
c3_hyp_silver_int_loc = find(diff2 == min(diff2))
Silver_branch_cut_loc_c3 = c3_imag(c3_hyp_silver_int_loc)

%% Plot Figure
figure('Name','Decay of the Creeping Wave');
plot(kxx,hyp_silver);
hold on
plot(c1_real, c1_imag);
plot(c3_real, c3_imag);
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

xlim([-10 10])
ylim([-10 10])

save nevels_michalski_creeping_wave_intersection.mat % Save data to plot the branch cuts

