%% Test Script
clearvars
% load('flap1_0Half_hydro.mat'); 
%% Constants
plotFlag = 0; 
rho      = 1000; 
g        = 9.81; 

tfOrder = 6; 
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
tfRadAR = cell(1,length(thetaN)); 
tfRadBR = cell(1,length(thetaN));
tfRadCR = cell(1,length(thetaN));
tfRadDR = cell(1,length(thetaN));

for i=1:length(thetaN)

filename = ["flap2_" + num2str(thetaN(i)) + "_hydro.mat"]; 
load(filename)

[~, A]      = getAddedMassFromBEMIO(hydro,plotFlag,rho);
[w,B]       = getRadiationDampingFromBEMIO(hydro,plotFlag,rho);
[tB,radIRF] = getRadiationIRFfromBEMIO(hydro,plotFlag,rho);
[~, E]      = getExcitationMagnitudeFromBEMIO(hydro,plotFlag,rho,g);
[tE,excIRF] = getExcitationIRFfromBEMIO(hydro,plotFlag,rho,g);

A_Pitch(:,i)    = A.Pitch;  % Added Mass
A_PitchInf(i)   = hydro.Ainf(5,5) * rho;  % Added Mass at infinite freq

KH_Pitch(i)     = hydro.Khs(5,5) * rho * g; % Hydrstatic Restoring Coefficient

B_Pitch(:,i)    = B.Pitch;  % Radiation Damping Coefficient
B_IRF(:,i)      = radIRF.Pitch;  % Radiation Damping Coefficient
B_IRFt(:,i)     = tB;  % Radiation Damping Coefficient
B_IRFSS(:,i)    = radIRF.PitchSS;  % Radiation Damping Coefficient

E_Pitch(:,i)    = E.Pitch;  % Excitation Force Magnitude
E_IRF(:,i)      = excIRF.Pitch;  % Excitation Damping Coefficient
E_IRFt(:,i)     = tE;  % Excitation Damping Coefficient

%% Getting the transfer function from the frequency response >>invfreqs 
% This is done by passing the freq. response and the angular 
% frequency vector through >>invfreqs 

% [bR(:,i), aR(:,i)] = invfreqs(B_Pitch(:,i),w,tfOrder,tfOrder);
% [bE(:,i), aE(:,i)] = invfreqs(E_Pitch(:,i),w,tfOrder,tfOrder);

% [bIRF_fft(:,i),freq] = Amp_Spectrum(B_IRF(:,i),1/(w(2)/2*pi),'rectangular','ECF');

%% Getting the transfer function from the frequency response >>tfest
% This is done by passing the freq. response and the angular 
% frequency vector through >>tfest
% sysR(i) = tfest(B_Pitch(:,i),tfOrder,tfOrder);
% sysE(i) = tfest(E_Pitch(:,i),tfOrder,tfOrder);

%% Getting the transfer function from the frequency response >>imp2ss
tSIRF = B_IRFt(2,i); 
[tfRadAR{i},tfRadBR{i},tfRadCR{i},tfRadDR{i},~,~] = imp2ss(B_IRF(:,i)*tSIRF,tSIRF);
% tfRad = ss2tf(AR3,BR3,CR3,DR3);
end

%% Save data

% new_filename = 'Flap2_DataALL_new'; 
% save(new_filename,'A_Pitch','A_PitchInf','B_Pitch','E_Pitch','KH_Pitch','w','B_IRF','B_IRFSS','B_IRFt','E_IRF','E_IRFt',...
%     'tfRadAR','tfRadBR','tfRadCR','tfRadDR')

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
