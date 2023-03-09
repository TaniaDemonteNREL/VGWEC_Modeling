function thetaDot = vgWecPlant(theta0, w_Wave, hydroCoeff, hydroVCoeff, Fexc, Frad)
% This function represents the plant of a variable geometry wave energy
% converter. The model here is derived by Mcabe 2006. 

% Inputs: 
% theta0      = initial state conditions
% w           = wave excitation angular frequency 'Omega' (rad/s)
% hydroCoeff  = time-invariable hydrodynamic coefficients
% hydrovCoeff = time-variable hydrodynamic coefficients
% Fexc        = excitation force from the waves (N)

% Outputs: 
% thetaDot(1)    = pitch angular velocity (rad/s). 
% thetaDot(2)    = pitch angular acceleration (rad/s^2). 

% States: 
theta1  = theta0(1); % (1) theta:  pitch angle response
theta2  = theta0(2); % (2) thetaD: pitch angle velocity

% Find the index for the discretized theta 
thetaN = find(deg2rad(hydroCoeff.theta) > abs(theta1),1); 

% Find the index for the wave frequency 
w = find(w_Wave > hydroCoeff.w,1); 

% thetaN = 1; 
% w  = 13; 

IQ = hydroCoeff.IQ(thetaN); % Moment of Inertia about the axis of pitch
IA = hydroVCoeff.IA(thetaN); % Added mass moment of inertia at inf. freq. (t)
% b  = hydroVCoeff.b(w, thetaN); % Radiation damping coefficient (t) - used
% in Frad instead
c  = hydroVCoeff.c(thetaN); % Hydrostatic restoring coefficient. (t) 

% Calculate Radiation Force
% Frad = fRad(theta0, w, hydroVCoeff, thetaN); 

% Spring Force to keep the flap from rotating
K = 5000; % Spring Constant
Fspring = K * theta1; 

% States Derivation
thetaDot     = zeros(1,2); 
thetaDot(1)  = theta2; 
thetaDot(2)  = (Fexc - Frad - c*theta1 - Fspring) / (IQ + IA); 

end 