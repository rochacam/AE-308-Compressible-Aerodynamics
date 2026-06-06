% Code Name: FFT Boundary Layer Aerodynamics Project
% Code Description: Reduces the hotwire data provided and use FFT to calculate turbulance kinetic energy and intentisy across boundary layer
% Author: Matheus Rocha Carlos
% Email: ROCHACAM@my.erau.edu
% Class: AE314 - Section 02DB
% Date: 04/06/2026 until 04/14/2026
% Worked With: N/a | Utilized Copilot1 (Ctrl+Shift+P) to generate plot codes and summary table
% FFT Code provided by: SignalAnalysisFourier_exported.pdf available in Canvas

%% LOAD DATA
% Iniciated clean MATLAB, imported all the data files, set sampling 
% frequency, constants and conversion factors for code. Extracted 
% calibration files, z-location file, and tare data required and sorted 
% them in respective array.
clear; clc; close all;
dataDir = fileparts(mfilename('fullpath')); 
fs = 20000; 
nu_air = 1.5e-5; 
R_air = 287; 
inHg2Pa = 3386.39; 
inWc2Pa = 249.089;
nWelch = 2^17; 
C = load(fullfile(dataDir,'Calibration.mat')); 
HW_cal = C.HW_cal(:)'; 
Pdyn_cal = C.Pdyn_cal(:)'; 
Patm_cal = C.Patm_cal(:)'; 
Tatm_cal = C.Tatm_cal(:)'; 
zMeas = C.zMeas_SI(:); 
Z = load(fullfile(dataDir, 'zero1.mat')); 
zData = Z.zero1;
V0_HW = mean(zData(:,1)); 
V0_Pdyn = mean(zData(:,2));

%% CALIBRATE ALL VOLATEGES TO RESPECTIVE UNITS USING CALIBARTION FILES
% Used calibration files and conversions to get Patm, Tatm, average density 
% and kinematic viscosity using sutherland law.
Patm_Pa = polyval(Patm_cal,mean(zData(:,3)))*inHg2Pa;
Tatm_K = polyval(Tatm_cal,mean(zData(:,4))); 
rho = Patm_Pa/(R_air*Tatm_K);
mu = 1.458e-6*Tatm_K^1.5/(Tatm_K+110.4);

%% IMPORT, SORT AND PROCESS EACH TRAVERSE POSITION FILES
% Found all data files and stored them into one array based on numbering of
% file. The imponted each extracted file into an empty array and using
% token for loop applied a new sort method to extract all total wall
% distances, full velocity versus time, and free-stream velocities. Then
% calculated mean and standard deviation of velocities and imported all
% data into array. And laslty calculate the free-stream velocity based on
% the mean of the free-stream velocity array.
listing = dir(fullfile(dataDir,'zPos_*.mat'));
fileNums = zeros(1,numel(listing)); 

% -------------------------------------------------------------------------
% For loop to tokenize each file number to make extraction easier, and sort 
% them by given data number.
% -------------------------------------------------------------------------
for k = 1:numel(listing)
    tok = regexp(listing(k).name,'zPos_(\d+)\.mat','tokens');
    fileNums(k) = str2double(tok{1}{1});
end
[fileNums,sortIdx] = sort(fileNums); 
listing = listing(sortIdx); 
nFiles  = numel(listing);
z_mm = zMeas(fileNums)*1e3;
U_mean = zeros(1,nFiles);
U_std = zeros(1,nFiles); 
U_inf_arr = zeros(1,nFiles);
U_all = cell(1,nFiles); 

% -------------------------------------------------------------------------
% for each token file extract the data as a matrix of 4 columns, tare the
% hot-wire voltage and dynamic preassure, and use calibration coefficients
% to get dynamic preassure, mean dynamic preassure, free-stream velocity,
% mean velocity, and standard deviation of velocity and store the data into
%  each respective pre-allocated array.
% -------------------------------------------------------------------------
for k = 1:nFiles 
    fname = fullfile(dataDir, listing(k).name); 
    D = load(fname);
    data = D.data;
    V_HW_tared = data(:,1)-V0_HW; 
    V_Pdyn_tared = data(:,2)-V0_Pdyn;
    U_hw = polyval(HW_cal,V_HW_tared);
    q_dyn = polyval(Pdyn_cal,V_Pdyn_tared)*inWc2Pa; 
    q_mean = mean(q_dyn);
    U_inf_k = sqrt(max(2*q_mean/rho,0));
    U_mean(k) = mean(U_hw); 
    U_std(k) = std(U_hw); 
    U_inf_arr(k) = U_inf_k; 
    U_all{k} = U_hw; 
end
U_inf = median(U_inf_arr);

%% CALCULATE BOUNDARY LAYER THICKNESS (Delta 99)
% Utilized linear interpolation to find where z has a local velocity to 
% free-stream ratio of 0.99, by creating an arrau of ratios of local 
% velocities to free-stream velocity. Then by sorting distances skipping 
% duplicates, and then sorting the ratio array by this distances and using 
% a find to get where sorted ratios equal to about 0.99 BL thickness we can 
% create the upward bracketing and downward bracketing points for 
% interpolation to get delta99.
u_ratio = abs(U_mean/U_inf); 
[z_sorted,si] = sort(z_mm); 
u_sorted = u_ratio(flip(si)); 
idx99 = find(u_sorted>=0.99,1,'first');
z1 = z_sorted(idx99-1); u1 = u_sorted(idx99-1); 
z2 = z_sorted(idx99); u2 = u_sorted(idx99);
delta99 = z1+(0.99-u1)*(z2-z1)/(u2-u1); 

%% CALCULATE WALL SHEAR STRESS AND SKIN FRICTION COEFFICIENT ASSUMING NO SLIP CONDITION
% Used only the most near-wall point, converting to mm, and 
% forcing local velocity of wall equal to 0, and get the velocity gradient
% through this points. Using the gradient calculate the shear stress, 
% friction velocity and skin friction coefficient.
nWall = min(1,nFiles);
z_wall = z_mm(1:nWall)'*1e-3; 
U_wall = [0;U_mean(1:nWall)'];
dUdz_wall = sum(z_wall.*U_wall)/sum(z_wall.^2);
tau_w = mu*dUdz_wall;
u_tau = sqrt(tau_w/rho);
Cf = tau_w/(0.5*rho*U_inf^2);

%% COMPARE FOUND BL DATA TO THE EMPIRICAL BOUNDARY LAYER PROFILES FOR COMPARISON
% Calculate eta(z/delta) and find velocity for a turbulent BL based on 
% 1/7 power law and use newtoniam method for Blausius laminar ODE to get 
% BL layout for laminar BL setting eta to known value for 99% velocity, 
% and then normalizing all etas to 0-1 range.
eta_fine = linspace(0, 1, 300);
U_turb = eta_fine.^(1/7); 
eta_B99 = 4.91; 
blasius_ode = @(xi,y)[y(2);y(3);-0.5*y(1)*y(3)];
opts = odeset('RelTol',1e-8,'AbsTol',1e-10); 
xispan = [0 eta_B99]; 
[xi_bl,Y_bl] = ode45(blasius_ode,xispan,[0 0 0.332],opts); 
U_lam = Y_bl(:,2);
eta_lam = xi_bl/eta_B99;

% -------------------------------------------------------------------------
% Plot the experimental boundary layer with error bars for the velocity at 
% each point wall height correcting all mean velocoties to postive values 
% for hotwire data and
% -------------------------------------------------------------------------
New_Umean = abs(U_mean);
fig1 = figure('Name','Fig 1 – Dimensional BL Profile','NumberTitle','off','Position',[50 50 650 520]);
errorbar(sortrows(New_Umean',"ascend"), z_mm, U_std', 'horizontal','r-', 'LineWidth',1.5, 'MarkerSize',6);
hold on;
xline(U_inf,'--c','LineWidth',1.2,'Label','U_{\infty}','LabelVerticalAlignment','bottom');
yline(delta99,'--g','LineWidth',1.2,'Label','delta99');
xlabel('Velocity  U  [m/s]',  'FontSize',13);
ylabel('Wall distance  z  [mm]', 'FontSize',13);
title('Dimensional Boundary Layer Profile', 'FontSize',14);
legend('Average velocity Profile', 'Location','northwest', 'FontSize',11);
grid on; box on; xlim([0, U_inf+5]); set(gca,'FontSize',12);

% -------------------------------------------------------------------------
% Plot the non-dimentsional Boundary layer and compare to laminar and 
% turbulent BL theorical data with no error bars, normalizing the measured 
% free-stream velocity mean to be positive and the eta up to delta99.
% -------------------------------------------------------------------------
fig2 = figure('Name','Fig 2 – Non-Dimensional BL Profile','NumberTitle','off','Position',[120 50 650 520]);
eta_meas = z_mm / delta99;
u_meas   = abs(U_mean / U_inf);
plot(sortrows(u_meas',"ascend"),eta_meas,'-b','MarkerSize',10,'LineWidth',2,'DisplayName','Experiment');
hold on;
plot(U_turb,eta_fine,'--r','LineWidth',2,'DisplayName','Turbulent (1/7 power law)');
plot(U_lam,eta_lam,'-.g','LineWidth',2,'DisplayName','Laminar (Blasius)');
xline(1,':y','LineWidth',1,'Label','z = delta99');
xlabel('\eta = z / delta99','FontSize',13);
ylabel('u / U_{\infty}','FontSize',13);
title('Non-Dimensional Boundary Layer Profile','FontSize',14);
legend('Location','southeast','FontSize',11);
grid on; box on; set(gca,'FontSize',12);
xlim([-0.1, max(eta_meas)*1.1]);  ylim([0, 1.4]);
% -------------------------------------------------------------------------
% From figure 2 it can be seen that the experimental boundary layer data
% overlays very well with the theorical turbulent boundary layer profile.
% This is indicated by the non-dimensaional velocity ratio rapid increase
% and delayed gradient with eta. This is expected as we also want to 
% calculate the turbulance intensity and PSD of the boundary layer 
% indicating we are supposed to analyze a tuerbulent boundary layer. 
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Plot the PSD Turbulent intensity and Tubulent Kinetic Energy using Welch.
% First selected the last measurement point that is the closest to 
% free-stream the most turbulent fluctuations and the near-wall for 
% position and comparison assuming a 50% overlap based on assumed welch.
% -------------------------------------------------------------------------
idxFS = nFiles; 
idxNW = 1; 
nOvlp = round(nWelch/2); 
fig3 = figure('Name','Fig 3 – Power Spectral Density via Welch Method');
locs = {idxFS, idxNW};
names = {'Outer region (z = %0.1f mm)', 'Near wall (z = %0.1f mm)'};

% -------------------------------------------------------------------------
% Build power spectral density then scale it to turbulent intensity 
% calculation and TKE. Done by getting locations of data, each 
% instantaneous velocities, the velocity fluctuations u' from the mean 
% velocity of each data, and use pwelch of u' to to aclculate u'_rms the 
% TU intensity and TKE assuming 1-component.
% -------------------------------------------------------------------------
for p = 1:2
    idx = locs{p};
    u_sig = U_all{idx}; 
    u_fl = u_sig-mean(u_sig);
    [Pxx,freq] = pwelch(u_fl,hann(nWelch),nOvlp,nWelch,fs);
    Tu_rms = sqrt(trapz(freq, Pxx)); 
    Tu = Tu_rms/U_mean(idx);
    TKE_spec = 0.5*Pxx;
    subplot(1,2,p); 
    loglog(freq, TKE_spec, '-b','LineWidth',1.2);
    ylim([10^(-10), 1]); xlim("tight"); hold on;

    % ---------------------------------------------------------------------
    % create Kolmogorov -5/3 reference line anchored at mid-range 
    % assumption and plot on log-log scale together with the welch 
    % experimental data
    % ---------------------------------------------------------------------
    f_ref = 500;
    idx_f = find(freq>=f_ref,1);
    K_amp = TKE_spec(idx_f)*f_ref^(5/3); 
    f_k = logspace(log10(1),log10(8000),100); 
    loglog(f_k,K_amp*f_k.^(-5/3),'--r','LineWidth',2,'DisplayName','-5/3 slope (Kolmogorov)');
    xlabel('Frequency [Hz]','FontSize',12);
    ylabel('TKE spectrum [m²/s²/Hz]','FontSize',12);
    title(sprintf(names{p}, z_mm(idx)),'FontSize',12);
    legend('TKE spectrum','-5/3 slope','Location','southwest','FontSize',10);
    grid on; box on; set(gca,'FontSize',11);
    text(0.05,0.92, sprintf('AVG Tu = %.1f%%', Tu*100),'Units','normalized','FontSize',11,'FontWeight','bold');
end
sgtitle('Power Spectral Density via Welch Method', 'FontSize',14);

% -------------------------------------------------------------------------
% Plot the turbulance intensity profile in the BL profile up to delta99 
% findiing turbulent porifle at each z from the standard deviation of u
% -------------------------------------------------------------------------
Tu_profile = U_std./U_mean;
fig4 = figure('Name','Fig 4 – Turbulent Intensity Profile','NumberTitle','off','Position',[330 50 600 480]);
plot(Tu_profile*100,z_mm,'r-','LineWidth',2,'MarkerSize',10);
yline(delta99,'--y','LineWidth',1.2,'Label','delta99');
xlabel('Turbulent Intensity Tu [%]','FontSize',13);
ylabel('Wall distance z [mm]','FontSize',13);
legend('Experimental Turbulance Intensity in BL','Location','northwest','FontSize',10);
title('Turbulent Intensity Profile', 'FontSize',14);
grid on; box on; set(gca,'FontSize',12);

%% SUMMARY TABLE OF ALL FINDINGS OF EXPERIMENTAL BL
fprintf('\n══════════════════════════════════════════════════\n');
fprintf('  BOUNDARY LAYER SUMMARY\n');
fprintf('══════════════════════════════════════════════════\n');
fprintf('  U_inf              = %8.3f m/s\n',U_inf);
fprintf('  Patm               = %8.4e \n',Patm_Pa);
fprintf('  Tatm               = %8.4e \n',Tatm_K);
fprintf('  rho                = %8.4f kg/m³\n', rho);
fprintf('  mu                 = %8.4e Pa·s\n',mu);
fprintf('  nu                 = %8.4e m²/s\n',mu/rho);
fprintf('  Delta99            = %8.3f mm\n',delta99);
fprintf('  dU/dz|_wall        = %8.1f s⁻¹\n',abs(mean(dUdz_wall)));
fprintf('  τ_w                = %8.4f Pa\n',abs(mean(tau_w)));
fprintf('  Cf                 = %8.4e \n',abs(mean(Cf)));
fprintf('══════════════════════════════════════════════════\n\n');