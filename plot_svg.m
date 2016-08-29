clear all;

%% lets plot 3 cycles of 50Hz AC voltage
f = 50;  % frequency
Vm = 10; % peak
phi = 0; % phase

% generate the signal
t = [0:0.0001:3/f];
th = 2*pi*f*t;
v = Vm*sin(th+phi);

% plot it
figure;
plot(t*1E3, v);

% change settings
opt = [];
opt.XLabel = 'Time, t (ms)';   % xlabel
opt.YLabel = 'Voltage, V (V)'; % ylabel
opt.Title = 'Voltage as a function of time'; % title
setPlotProp(opt);
opt.Colors = [0, 0, 0];   % [red, green, blue]
opt.LineWidth = 2;        % line width
opt.LineStyle = {'--'};   % line style: '-', ':', '--' etc