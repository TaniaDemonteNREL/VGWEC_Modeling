%% VG WEC Model Setup
% This is the second version where I am incorporating real data. 
%% Clean up & Setup

clearvars; close all; clc;
load('Flap1_Data.mat');

%% Sim Setup

dt   = 0.01; % Simulation time step (s) 
tSim = 100; % Simulation time (s)
t    = 0:dt:tSim; % Time vector (s)

dTheta = 10; % Discretized Theta Increments in degrees

% Initial Conditions
theta0 = [0,0];  

%% Temp Info

I55 = 1.3754e+05 * 0.78; % Pitch mass moment of inertia Tom. N 2016 Eq(1)
C55 = 2.0233e+05 * 0.78; % Hydrostatic Restoring Coefficient Tom. N 2016 Eq(2)

%% Time-Invariable Hydrodynamic Coefficients
hydroCoeff.w     = w; % wave excitation angular frequency (rad/s)
hydroCoeff.IQ    = I55 * ones(length(A_Pitch),1); % Moment of Inertia about the axis of pitch
hydroCoeff.theta = 0:dTheta:90;  
thetaN           = length(hydroCoeff.theta); % Number of discretized points. Every 10 degrees.
%% Time-Variable Hydrodynamic Coefficients

hydroVCoeff.IA  =  A_Pitch; % Added mass moment of inertia at inf. freq. (t)
hydroVCoeff.b   =  B_Pitch; % Radiation damping coefficient (t) 
% hydroVCoeff.c   =  KH_Pitch; % Hydrostatic restoring coefficient. (t) 
hydroVCoeff.c   =  randi([round(C55*0.9) round(C55*1.1)], length(hydroCoeff.w),thetaN); % Hydrostatic restoring coefficient. (t) 
%% External Forces

rho_g_h2_w = 2452500; 
hydroVCoeff.KE =  E_Pitch;
hydroVCoeff.KR =  B_Pitch;

%% Wave Excitation

w_Wave   = 0.5; % Wave excitation angular frequency 'Omega' (rad/s)
amp_Wave = 0.2; % Wave excitation amplitude (m)
