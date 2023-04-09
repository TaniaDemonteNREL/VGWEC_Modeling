function FE = fExcInterp(theta0, eta, w_Wave,hydroVCoeff,hydroCoeff)
% This function generates the linear excitation force

% Inputs: 
% theta0      = current states condition
% eta         = wave excitation (m)
% w_Wave      = wave excitation angular frequency 'Omega' (rad/s)
% hydroCoeff  = time-invariable hydrodynamic coefficients
% hydrovCoeff = time-variable hydrodynamic coefficients


% Outputs: 
% FE    = excitation force (N) 
theta1 = theta0(1); 

% Find the index for the discretized theta 
theta1N = wrapToPi(theta1); 

if abs(rad2deg(theta1N)) >= 95 
    theta1N = pi - abs(theta1N); 
end 

theta1Ndeg = rad2deg(theta1N);

thetaN = find(deg2rad(hydroCoeff.fineTheta) > abs(theta1N),1)-1; 


% Find the index for the wave excitation frequency 
w  = find(hydroCoeff.w > w_Wave,1)-1;


KE = hydroVCoeff.FExcInterp(w, thetaN); % Added mass moment of inertia at inf. freq. (t)

FE = zeros(1,1); 
FE(1) = eta*KE; 