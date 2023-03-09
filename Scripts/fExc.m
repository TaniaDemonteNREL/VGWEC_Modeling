function FE = fExc(theta0, eta, w_Wave,hydroVCoeff,hydroCoeff)
% This function generates the linear excitation force

% Inputs: 
% theta0      = current states condition
% eta         = wave excitation (m)
% w_Wave      = wave excitation angular frequency 'Omega' (rad/s)
% hydroCoeff  = time-invariable hydrodynamic coefficients
% hydrovCoeff = time-variable hydrodynamic coefficients


% Outputs: 
% FE    = excitation force (N) 

% Find the index for the discretized theta 
thetaN = find(deg2rad(hydroCoeff.theta) > abs(theta0(1)),1)-1; 

% Find the index for the wave excitation frequency 
w  = find(hydroCoeff.w > w_Wave,1);

% thetaN = 1; 
% w  = 13; 

KE = hydroVCoeff.KE(w, thetaN); % Added mass moment of inertia at inf. freq. (t)

FE = zeros(1,1); 
FE(1) = eta*KE; 
