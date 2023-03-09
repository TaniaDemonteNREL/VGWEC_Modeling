function FR = fRadSim(theta0, w_Wave,hydroVCoeff,hydroCoeff)
% This function generates the radiation force

% Inputs: 
% theta0      = current states condition
% w           = wave frequency index
% hydrovCoeff = time-variable hydrodynamic coefficients
% thetaN      = pitch angle index

% Outputs: 
% FR    = radiation force (N) 

% Find the index for the discretized theta 
thetaN = find(deg2rad(hydroCoeff.theta) > abs(theta0(1)),1)-1; 

% Find the index for the wave excitation frequency 
w  = find(hydroCoeff.w > w_Wave,1);

% thetaN = 1; 
% w  = 13; 

theta2  = theta0(2); 
b       = hydroVCoeff.b(w, thetaN); % Radiation damping coefficient (t) 
KR      = b; % For single harmonic wave

FR    = zeros(1,1); 
FR(1) = theta2*KR; 