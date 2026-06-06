% Code Name: Heat Transfer rate calculator
% Code Description: Axisymetric missile noseheat transfer rate fro wall
% temperatures
% Author: Matheus Rocha Carlos
% Email: ROCHACAM@my.erau.edu
% Class: AE308 - Section 1
% Date: 11/15/2025
% Worked With: N/a

% Iniciation
clear; clc;

% Given data and constant
gamma = 1.28;
R = 287;                 
T_inf = 288.15;           
rho_inf = 1.225;        
M = 8;
R_n= 0.005;        % in meters
eps = 1.0;                % emissivity
sigma = 5.670374419e-8;   % Stefan-Boltzmann constant
c1 = 1.83e-4;             % Sutton-Graves constant in meters

% Flow properties calculation
a_inf = sqrt(gamma*R*T_inf);
V = M*a_inf;

T0 = T_inf*(1 + (gamma-1)/2*M^2);      % stagnation temperature
cp = gamma/(gamma-1)*R;
h0 = cp*T0;                             % total enthalpy J/kg

fprintf("Stagnation temperature T0 = %.2f K\n", T0);
fprintf("Total enthalpy h0 = %.3e J/kg (%.3f MJ/kg)\n", h0, h0/1e6);

% Convective heating calculation
q_conv = (c1/sqrt(R_n) * (1 - ((cp*T_inf)/h0))) * sqrt(rho_inf) * V^3;   % W/m^2

fprintf("Convective heat flux q_conv = %.3e W/m^2 (%.2f MW/m^2)\n", ...
         q_conv, q_conv/1e6);

% Radiative cooling calculation
Twall = linspace(300,1600,30000);
q_rad = eps * sigma * Twall.^4;
q_net = q_conv - q_rad;

% Plot results of nose heat transfer
figure; hold on; grid on;
plot(Twall, q_conv*ones(size(Twall)), 'LineWidth', 2);
plot(Twall, q_rad, 'LineWidth', 2);
plot(Twall, q_net, 'LineWidth', 2);
xline(1600, '--k', '1600 K = Vaporization', 'LabelVerticalAlignment','bottom');
xlabel('Wall Temperature T_w (K)');
ylabel('Heat Flux (W/m^2)');
legend('Convective (SG)', 'Radiative (εσT^4)', 'Net Heat Flux',Location='southwest');
title('Stagnation-Point Heating vs Wall Temperature (1-cm Nose, M = 8)');


%% Tungsten Nose Lumped Transient Heating to find time to vaporization
clear; clc;

% Given
rho = 19250;          % kg/m^3
cp = 135;             % J/kg-K
L = 0.005;            % cube side in meters (0.5 cm)
A = 6*L^2;            % surface area
Vcube = L^3;          % cube volume
m = rho * Vcube;      % mass of tungsten nose

% Convective heat flux (Chapman/Sutton–Graves)
gamma = 1.28;
R = 287; T_inf = 288.15; rho_inf = 1.225; M = 8;
c1 = 1.83e-4;
eps = 1.0;
sigma = 5.670374419e-8;

a_inf = sqrt(gamma*R*T_inf);
V = M*a_inf;
q_conv = (c1/sqrt(L) * (1 - ((cp*T_inf)/(cp*T_inf*(V^2/2))))) * sqrt(rho_inf) * V^3;  % W/m^2

% Time stepping
dt = 0.0001;                 % s
T = 300;                   % initial temp
T_target = 1600;           % K
time = 0;

T_history = T;
t_history = time;

while T < T_target
    q_rad = eps * sigma * T^4;
    q_net = (q_conv - q_rad);       % W/m^2
    Q_net = q_net * A;              % total W
    
    dT = Q_net*dt / (m*cp);         % K
    T = T + dT;
    
    time = time + dt;
    T_history(end+1) = T;
    t_history(end+1) = time;
end

fprintf('Time to reach %.0f K: %.4f seconds\n', T_target, time);

% Plot temperature vs time
figure;
plot(t_history, T_history, 'LineWidth',2);
xlabel('Time (s)');
ylabel('Tungsten Nose Temperature (K)');
grid on;
title('Transient Lumped-Element Heating of Tungsten Nose');

