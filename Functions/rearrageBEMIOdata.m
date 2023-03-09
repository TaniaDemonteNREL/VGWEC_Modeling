%% Test Script

% load('flap1_0Half_hydro.mat'); 
%% Constants
plotFlag = 0; 
rho      = 1000; 
g        = 9.81; 

%% Form vectors at each discretized theta

dTheta = 10; % Discretized Theta Increments in degrees

thetaN = 0:dTheta:90; 

%% Memory allocation
A_Pitch     = zeros(500,length(thetaN));
A_PitchInf  = zeros(length(thetaN),1);
B_Pitch  = zeros(500,length(thetaN));
% radIRF = zeros(length(thetaN),1);
E_Pitch  = zeros(500,length(thetaN));
KH_Pitch = zeros(length(thetaN),1);

for i=1:length(thetaN)

filename = ["flap2_" + num2str(thetaN(i)) + "_hydro.mat"]; 
load(filename)

[~, A]      = getAddedMass(hydro,plotFlag,rho);
[w,B]       = getRadiationDamping(hydro,plotFlag,rho);
% [tt,radIRF] = getRadiationIRF(hydro,plotFlag,rho);
[~, E]      = getExcitationMagnitude(hydro,plotFlag,rho,g);

A_Pitch(:,i)    = A.Pitch;  % Added Mass
A_PitchInf(i)   = hydro.Ainf(5,5) * rho;  % Added Mass at infinite freq
B_Pitch(:,i) = B.Pitch;  % Radiation Damping Coefficient
E_Pitch(:,i) = E.Pitch;  % Excitation Force Magnitude
KH_Pitch(i)  = hydro.Khs(5,5) * rho * g; % Hydrstatic Restoring Coefficient

end

save('Flap2_Data','A_Pitch','A_PitchInf','B_Pitch','E_Pitch','KH_Pitch','w')

%%
% AA_Surge = squeeze(hydro.ss_A(1,1,:,:)); 
% AA_Heave = squeeze(hydro.ss_A(2,2,:,:));
% AA_Pitch = squeeze(hydro.ss_A(3,3,:,:));
% 
% BB_Surge = squeeze(hydro.ss_B(1,1,:));
% BB_Heave = squeeze(hydro.ss_B(2,2,:));
% BB_Pitch = squeeze(hydro.ss_B(3,3,:));
% 
% CC_Surge = squeeze(hydro.ss_C(1,1,:,:));
% CC_Heave = squeeze(hydro.ss_C(2,2,:,:));
% CC_Pitch = squeeze(hydro.ss_C(3,3,:,:));
% CC_Heave = CC_Heave';
% CC_Surge = CC_Surge';
% CC_Pitch = CC_Pitch';
% 
% DD_Heave = squeeze(hydro.ss_D(2,2));
% DD_Surge = squeeze(hydro.ss_D(1,1));
% DD_Pitch = squeeze(hydro.ss_D(3,3));
