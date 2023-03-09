function FR = fRad(theta0, w,hydroVCoeff,thetaN)
% This function generates the radiation force

% Inputs: 
% theta0      = current states condition
% w           = wave frequency index
% hydrovCoeff = time-variable hydrodynamic coefficients
% thetaN      = pitch angle index

% Outputs: 
% FR    = radiation force (N) 

theta2  = theta0(2); 
KR      = hydroVCoeff.KR(w, thetaN); % Added mass moment of inertia at inf. freq.

FR    = zeros(1,1); 
FR(1) = theta2*KR; 