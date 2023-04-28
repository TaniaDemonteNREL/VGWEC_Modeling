function [AR, BR, CR, DR] = getRadTFCoeffs(theta0,tfRadAR_ROm, tfRadBR_ROm, tfRadCR_ROm, tfRadDR_ROm, hydroCoeff)
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

% AR       = tfRadAR{thetaN}; % A Matrix Coefficients
% BR       = tfRadBR{thetaN}; % A Matrix Coefficients
% CR       = tfRadCR{thetaN}; % A Matrix Coefficients
% DR       = tfRadDR{thetaN}; % A Matrix Coefficients

% AR       = tfRadAR_RO{thetaN}; % A Matrix Coefficients
% BR       = tfRadBR_RO{thetaN}; % A Matrix Coefficients
% CR       = tfRadCR_RO{thetaN}; % A Matrix Coefficients
% DR       = tfRadDR_RO{thetaN}; % A Matrix Coefficients

% Get transfer function coefficients depending on the discretized angle. 
AR       = zeros(6,6); 
BR       = zeros(6,1); 
CR       = zeros(1,6); 
DR       = 0; 

AR       = tfRadAR_ROm(:,(thetaN-1)*6+1:(thetaN*6)); % A Matrix Coefficients
BR       = tfRadBR_ROm(:,thetaN); % A Matrix Coefficients
CR       = tfRadCR_ROm(1,(thetaN-1)*6+1:(thetaN*6)); % A Matrix Coefficients
DR       = tfRadDR_ROm(thetaN); % A Matrix Coefficients