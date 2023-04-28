% Support Script 

d = 10; % Water Depth (m) 

h = 10;   % OSWEC Height (m)
t = 0.33; % OSWEC Thickness (m)
w = 5;    % OSWEC Width 
V = 14.1; % OSWEC Volume (m^3)

tf = 0.33; % Flap minor axis (m)
hf = 2;    % Flap major axis (m)
wf = 4.5;  % Flap width (m)
ws = 0.25; % Side support width (m)
rb = 4.8; % Radial Center of Buoyancy (m)

rhoW = 1000; % Density of water (kg/m^3) 
rhoD = rhoW/2; % Density of the device (kg/m^3)  
g    = 9.81; % Gravitational Acceleration (m/s^2)

%% Calculations

% Pitch mass moment of inertia Eq(1)
I55 = rhoD*w*t*h^3/6 *(1+(t/(2*h))^2)

% Hydrostatic Restoring Coefficient Eq(2)
c55 = rhoD*w*t*h^2*g/4
