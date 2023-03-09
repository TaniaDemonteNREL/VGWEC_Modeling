function [bR, aR] = getRadTFCoeffs(theta0,hydroVCoeff,hydroCoeff)
% This function extracts the coefficients of the Radiation Transfer
% Function based on previously fitted data between discretized thetas. 

% Inputs: 
% theta0      = current states condition
% hydrovCoeff = time-variable hydrodynamic coefficients
% thetaN      = pitch angle index

% Outputs: 
% bR    = radiation force (N) 

% Find the index for the discretized theta 
thetaN = find(deg2rad(hydroCoeff.theta) > abs(theta0(1)),1)-1; 

% Get transfer function coefficients depending on the discretized angle. 

bR       = hydroVCoeff.bR(thetaN); % Numerator Coefficients
aR       = hydroVCoeff.aR(thetaN); % Denominator Coefficients
