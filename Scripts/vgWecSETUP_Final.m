%% VG WEC Model Setup
% This is the main file to run the time-varying simulation. Simply run this
% file to get started. 
% The TV model can be selected by changing the Time_Variant variable in
% line 23, and the regular wave frequency and amplitude can be changed in
% line 26 and 29 respectevely. 

% Author: Tania Demonte Gonzalez email: tsdemont@mtu.edu
% Date:   4/28/2023
%% Clean up & Setup

clearvars; close all; clc;
load('Flap2_DataAll.mat'); % WAMIT Data
load('B1Data.mat'); % 
load('flap2_Frad_Reduced.mat')
rng('default')
%% Time - Varying Subsystem 
% Time_Variant == 0 for a LTI hydrodynamics model 
% Time_Variant == 1 for the LTV hydrodynamics model - Discrete 
% Time_Variant == 2 for the LTV hydrodynamics model - Interpolated 
% Time_Variant == 3 for the LTV hydrodynamics model - Polyfit 

Time_Variant = 0;  
%% Wave Excitation Parameters

w_Wave   = 0.8; % Wave excitation angular frequency 'Omega' (rad/s)
T_wave   = (2*pi)/w_Wave; % Wave Period (s)

amp_Wave = 0.5; % Wave excitation amplitude (m)

w_WaveIndx           = find(w > w_Wave,1); 
hydroCoeff.Kspring   = 10000; % Spring coefficient

%% Sim Setup

dt    = 0.01; % Simulation time step (s) 
tSim  = 300; % Simulation time (s)
t     = 0:dt:tSim; % Time vector (s)
tRamp = 100; 

dTheta = 10; % Discretized Theta Increments in degrees / Only change if you have new WAMIT data to add

% Initial Conditions
theta0 = [0,0];  

%% Time-Invariable Hydrodynamic Coefficients
hydroCoeff.w     = w; % wave excitation angular frequency (rad/s)
hydroCoeff.IQ    = 672.72; % Moment of Inertia about the axis of pitch from Solidworks Px
hydroCoeff.theta = 5:dTheta:95;  
thetaN           = length(hydroCoeff.theta); % Number of discretized points. Every 10 degrees.

hydroCoeff.singleThetaIndex = 1; % Coefficient corresponding to a fully-closed flap for the LTI model
%% Time-Variable Hydrodynamic Coefficients

hydroVCoeff.IA  =  A_PitchInf; % Added mass moment of inertia at inf. freq. (t)
hydroVCoeff.b   =  B_Pitch; % Radiation damping coefficient (t) 
hydroVCoeff.c   =  KH_Pitch; % Hydrostatic restoring coefficient. (t) 

%% External Forces

hydroVCoeff.KE =  E_Pitch;
hydroVCoeff.KR =  B1_all;

%% Format Rad Data 
% This data is obtained using imp2ss in rearrangeBemio.m and their order is
% reduced in ModelReduction_RadiationSS.m

tfRadAR_ROm = cell2mat(tfRadAR_RO); 
tfRadBR_ROm = cell2mat(tfRadBR_RO); 
tfRadCR_ROm = cell2mat(tfRadCR_RO); 
tfRadDR_ROm = cell2mat(tfRadDR_RO); 

%% Run Interpolation and Polyfit Script

new_dTheta = 1; % Delta Theta to be new_dTheta degree increments for the interpolation
polyOrder  = 6; % Order of the polyfit function
run('Scripts\dataInterpolation.m')

%% Radiation SS IMP2SS METHOD 

% tSIRF       = B_IRFt(2,hydroCoeff.singleThetaIndex); 
% sys3        = imp2ss(B_IRF(:,hydroCoeff.singleThetaIndex)*tSIRF,tSIRF);
% 
% ssOrder     = 6; 
% sys_reduced = balred(sys3,ssOrder);
% 
% AR3 = sys_reduced.A;  
% BR3 = sys_reduced.B; 
% CR3 = sys_reduced.C; 
% DR3 = sys_reduced.D; 

%% New rad SS from SS fitting toolbox
% load('C:\Users\tdemonte\OneDrive - NREL\VGWEC\Full Flap Data\Flap2_FDI_Frad.mat')
% 
% AR3 = radA{hydroCoeff.singleThetaIndex};  
% BR3 = radB{hydroCoeff.singleThetaIndex}; 
% CR3 = radC{hydroCoeff.singleThetaIndex}; 
% DR3 = zeros(1,1);

%% Plot Simulation Output

% simOut = sim('VGWEC_Model.slx'); 
% 
% %% Get Sim Data 
% 
% sim.fExc     = simOut.logsout{1}.Values; 
% sim.eta      = simOut.logsout{2}.Values;  
% sim.theta    = simOut.logsout{3}.Values; 
% sim.fRad     = simOut.logsout{4}.Values;
% sim.thetaDot = simOut.logsout{5}.Values; 
% sim.thetaDeg = simOut.logsout{6}.Values;   
% 
% %% Plot
% LW = 'LineWidth'; 
% 
% figure; subplot(3,1,1); plot(sim.thetaDeg,LW,1.5); hold on; grid on
% xlabel('Time (s)'); ylabel('Flap Response');
% legend('Angle Displacement (degrees)' , 'Angular Velocity (degrees/s)');
% title('Pitch Flap Response')
% 
% subplot(3,1,2); plot(sim.fExc,'b',LW,1.5); grid on
% xlabel('Time (s)'); ylabel('Excitation Force (N)')
% title('Excitation Force')
% subplot(3,1,3); plot(sim.fRad,'r',LW,1.5); grid on
% xlabel('Time (s)'); ylabel('Radiation Force (N)')
% title('Radiation Force')