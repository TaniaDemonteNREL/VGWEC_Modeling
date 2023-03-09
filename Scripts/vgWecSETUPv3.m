%% VG WEC Model Setup
%% Clean up & Setup

clearvars; close all; clc;
Simulink.sdi.clear; 
load('Flap2_Data.mat'); % WAMIT Data
load('Flap2_Ixx.mat'); % SolidWorks Data
rng('default')

%% Sim Setup

dt   = 0.01; % Simulation time step (s) 
tSim = 100; % Simulation time (s)
t    = 0:dt:tSim; % Time vector (s)

dTheta = 10; % Discretized Theta Increments in degrees

% Initial Conditions
theta0 = [0,0];  

%% Temp Info

I55 = 8726.88; % From Ixx in Solidworks
C55 = 2.0233e+05 * 0.78; % Hydrostatic Restoring Coefficient Tom. N 2016 Eq(2)

%% Time-Invariable Hydrodynamic Coefficients
hydroCoeff.w     = w; % wave excitation angular frequency (rad/s)
hydroCoeff.IQ    = Ixx(:,2); % Moment of Inertia about the axis of pitch
hydroCoeff.theta = 0:dTheta:90;  
thetaN           = length(hydroCoeff.theta); % Number of discretized points. Every 5 degrees.
%% Time-Variable Hydrodynamic Coefficients

hydroVCoeff.IA  =  A_PitchInf; % Added mass moment of inertia at inf. freq. (t)
hydroVCoeff.b   =  B_Pitch; % Radiation damping coefficient (t) 
hydroVCoeff.c   =  KH_Pitch; % Hydrostatic restoring coefficient. (t) 
% hydroVCoeff.c   =  randi([round(C55*0.95) round(C55*1.05)], length(hydroCoeff.w),thetaN); % Hydrostatic restoring coefficient. (t) 

%% External Forces

hydroVCoeff.KE =  E_Pitch;
hydroVCoeff.KR =  B_Pitch;

%% Wave Excitation

w_Wave   = 0.5; % Wave excitation angular frequency 'Omega' (rad/s)
amp_Wave = 1; % Wave excitation amplitude (m)

%% Plot Simulation Output

simOut = sim('VGWEC_Model.slx'); 

%% Get Sim Data 

sim.fExc     = simOut.logsout{1}.Values; 
sim.eta      = simOut.logsout{2}.Values;  
sim.theta    = simOut.logsout{3}.Values; 
sim.fRad     = simOut.logsout{4}.Values;
sim.thetaDot = simOut.logsout{5}.Values; 
sim.thetaDeg = simOut.logsout{6}.Values;   

%% Plot
LW = 'LineWidth'; 

figure; subplot(3,1,1); plot(sim.thetaDeg,LW,1.5); hold on; grid on
xlabel('Time (s)'); ylabel('Flap Response');
legend('Angle Displacement (degrees)' , 'Angular Velocity (degrees/s)');
title('Pitch Flap Response')

subplot(3,1,2); plot(sim.fExc,'b',LW,1.5); grid on
xlabel('Time (s)'); ylabel('Excitation Force (N)')
title('Excitation Force')
subplot(3,1,3); plot(sim.fRad,'r',LW,1.5); grid on
xlabel('Time (s)'); ylabel('Radiation Force (N)')
title('Radiation Force')