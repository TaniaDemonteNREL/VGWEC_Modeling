%% VG WEC Model Setup
%% Clean up & Setup

%clearvars; close all; clc;
rng('default')

%% Sim Setup

dt   = 0.01; % Simulation time step (s) 
tSim = 100; % Simulation time (s)
t    = 0:dt:tSim; % Time vector (s)

dTheta = 5; % Discretized Theta Increments in degrees

% Initial Conditions
theta0 = [0,0];  

%% Temp Info

I55 = 1.3754e+05 * 0.78; % Pitch mass moment of inertia Tom. N 2016 Eq(1)
C55 = 2.0233e+05 * 0.78; % Hydrostatic Restoring Coefficient Tom. N 2016 Eq(2)

%% Time-Invariable Hydrodynamic Coefficients
hydroCoeff.w     = 0.2:0.1:5; % wave excitation angular frequency (rad/s)
hydroCoeff.IQ    = randi([I55*0.9 I55*1.1], length(hydroCoeff.w),1); % Moment of Inertia about the axis of pitch
hydroCoeff.theta = 0:dTheta:90;  
thetaN           = length(hydroCoeff.theta); % Number of discretized points. Every 5 degrees.
%% Time-Variable Hydrodynamic Coefficients

hydroVCoeff.IA  =  randi([I55*1 I55*10], length(hydroCoeff.w),thetaN); % Added mass moment of inertia at inf. freq. (t)
hydroVCoeff.b   =  randi([I55.*0.5 I55*3], length(hydroCoeff.w),thetaN); % Radiation damping coefficient (t) 
hydroVCoeff.c   =  randi([round(C55*0.9) round(C55*1.1)], length(hydroCoeff.w),thetaN); % Hydrostatic restoring coefficient. (t) 

%% External Forces

rho_g_h2_w = 2452500; 
hydroVCoeff.KE =  randi([round(rho_g_h2_w*0.03) round(rho_g_h2_w*0.07)], length(hydroCoeff.w),thetaN);
hydroVCoeff.KR =  randi([round(rho_g_h2_w*0.003) round(rho_g_h2_w*0.007)], length(hydroCoeff.w),thetaN);

%% Wave Excitation

w_Wave   = 0.5; % Wave excitation angular frequency 'Omega' (rad/s)
amp_Wave = 0.2; % Wave excitation amplitude (m)

%% Troubleshooting

%  eta   = amp_Wave * sin(w_Wave * t); 
%  theta00 = [10*pi/180, 15*pi/180]; 
% 
%  theta(1,:) =  theta00(1) * sin(w_Wave * t) + theta00(1); 
%  theta(2,:) =  theta00(2) * cos(w_Wave * t) + theta00(2); 
 %%
%  clear FE FR
% 
%  for i = 1:length(t) 
%  FE(i) = fExc(theta(:,i), eta(i), w_Wave,hydroVCoeff,hydroCoeff); 
%  FR(i) = fRad(theta(:,i), w_Wave,hydroVCoeff,hydroCoeff); 
% 
%  end
