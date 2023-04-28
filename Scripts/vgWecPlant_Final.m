function thetaDot = vgWecPlant_Final(theta0, hydroCoeff, hydroVCoeff, Fhydro,Fpto, Time_Variant)
% This function represents the plant of a variable geometry wave energy
% converter. The model here is derived by Mcabe 2006. 

% Inputs: 
% theta0      = initial state conditions
% w           = wave excitation angular frequency 'Omega' (rad/s)
% hydroCoeff  = time-invariable hydrodynamic coefficients
% hydrovCoeff = time-variable hydrodynamic coefficients
% Fexc        = excitation force from the waves (N)

% Time_Variant == 0 for a LTI hydrodynamics model 
% Time_Variant == 1 for the LTV hydrodynamics model - Discrete 
% Time_Variant == 2 for the LTV hydrodynamics model - Interpolated 
% Time_Variant == 3 for the LTV hydrodynamics model - Polyfit 


% Outputs: 
% thetaDot(1)    = pitch angular velocity (rad/s). 
% thetaDot(2)    = pitch angular acceleration (rad/s^2). 

% States: 
theta1  = theta0(1); % (1) theta:  pitch angle response
theta2  = theta0(2); % (2) thetaD: pitch angle velocity

% Find the index for the discretized theta 
theta1N = wrapToPi(theta1); 
    

if Time_Variant == 1 % Discrete LTV Model
    if abs(rad2deg(theta1N)) >= 95 
        theta1N = pi - abs(theta1N); 
    end 
    thetaN = find(deg2rad(hydroCoeff.theta) > abs(theta1N),1);    
    IA     = hydroVCoeff.IA(thetaN); % Added mass moment of inertia at inf. freq. (t)
    %disp('Discrete LTV Model')

elseif Time_Variant == 2 % Interpolated LTV Model
    if abs(rad2deg(theta1N)) >= 95 
        theta1N = pi - abs(theta1N); 
    end 
    thetaN = find(deg2rad(hydroCoeff.fineTheta) > abs(theta1N),1);    
    IA     = hydroVCoeff.AinfInterp(thetaN); % Added mass moment of inertia at inf. freq. (t)
    %disp('Interpolated LTV Model')

elseif Time_Variant == 3 % Polyfit LTV Model
    theta1D = rad2deg(theta1); 
    if abs(theta1D) >= 95
        theta1D = 180 - abs(theta1D);
    end 
    IA = polyval(hydroVCoeff.AinfFit, abs(theta1D),[], hydroVCoeff.a_mu); 
    %disp('Polyfit LTV Model')

else % LTI Model
    IA = hydroVCoeff.IA(hydroCoeff.singleThetaIndex); % Added mass moment of inertia at inf. freq. (t) 
    %disp('LTI Model')

end 

IQ     = hydroCoeff.IQ; % Moment of Inertia about the axis of pitch (Constant)

% States Derivation
thetaDot     = zeros(1,2); 
thetaDot(1)  = theta2; 
thetaDot(2)  = (Fhydro - Fpto) / (IQ + IA); % thetaDot(2)  = (Fexc - Frad - Fbuoy - Fspring) / (IQ + IA); 
end 