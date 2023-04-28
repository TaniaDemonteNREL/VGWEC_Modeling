function thetaDot = vgWecPlantInterp(theta0, hydroCoeff, hydroVCoeff, Fhydro,Fpto, Time_Variant)
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
theta1  = theta(1); % (1) theta:  pitch angle response
theta2  = theta(2); % (2) thetaD: pitch angle velocity


% Find the index for the discretized theta 
theta1N = wrapToPi(theta1); 

if abs(rad2deg(theta1N)) >= 90 
    theta1N = pi - abs(theta1N); 
end 

if Time_Variant > 0
    thetaN = find(deg2rad(hydroCoeff.fineTheta) > abs(theta1N),1)-1;    
    IA     = hydroVCoeff.AinfInterp(thetaN); % Added mass moment of inertia at inf. freq. (t)
else
    IA = hydroVCoeff.IA(hydroCoeff.singleThetaIndex); % Added mass moment of inertia at inf. freq. (t) 
end 

IQ     = hydroCoeff.IQ; % Moment of Inertia about the axis of pitch

% States Derivation
thetaDot     = zeros(1,2); 
thetaDot(1)  = theta2; 
thetaDot(2)  = (Fhydro - Fpto) / (IQ + IA); % thetaDot(2)  = (Fexc - Frad - Fbuoy - Fspring) / (IQ + IA); 
end 