%% Radiation Damping
% This script is to generate the convolution part of the radiation force coefficients. 
% Different methods will be used and compared. 
% clearvars
% load('Flap2_Data.mat') 
dTheta = 10; % Discretized Theta Increments in degrees
thetaN = 0:dTheta:90; 
tfOrder = 6; 
phi = zeros(length(B_Pitch),1); 

for i=1:length(thetaN)

%% Getting the transfer function from the frequency response
% This is done by passing the freq. response and the angular 
% frequency vector through >>invfreqs 

[bIRF_fft(:,i),freq(:,i)] = Amp_Spectrum(B_IRF(:,i),1/(2*pi/w(1)),'rectangular','ECF');

[bR(:,i), aR(:,i)] = invfreqs(bIRF_fft(:,i),freq(:,i),tfOrder,tfOrder);
% [bE(:,i), aE(:,i)] = invfreqs(E_Pitch(:,i).*exp(1j*phi),w,tfOrder,tfOrder);

% [bIRF_fft(:,i),freq] = Amp_Spectrum(B_IRF(:,i),1/(w(2)/2*pi),'rectangular','ECF');

% This is done by passing the freq. response and the angular 
% frequency vector through >>tfest
% sysR(i) = tfest(B_Pitch(:,i),tfOrder,tfOrder);
% sysE(i) = tfest(E_Pitch(:,i),tfOrder,tfOrder);

% [ss,totbnd(i),hsv(:,i)] = imp2ss(B_IRF);
end

%% Transfer Function and SS
thetaIndx = 10; 

KradTF1 = tf(bR(:,thetaIndx)',aR(:,thetaIndx)'); 
[AR1,BR1,CR1,DR1] = tf2ss(bR(:,thetaIndx)',aR(:,thetaIndx)');

% % Test Result
% [y_New1,t_New1]=impulse(KradTF1);
% 
% %
% figure; 
% plot(B_IRFt(:,thetaIndx),B_IRF(:,thetaIndx)); hold on;
% plot(t_New1,y_New1)

%% Toolbox TF 
load('KradTF_Toolbox.mat')
KradTF2 = Krad; 
[AR2,BR2,CR2,DR2] = tf2ss(KradTF2.Numerator{1,1},KradTF2.Denominator{1,1});
%% IMP2SS METHOD 
tSIRF = B_IRFt(2,10); 
[AR3,BR3,CR3,DR3,totbnd,hsv] = imp2ss(B_IRF(:,10)*tSIRF,tSIRF);
sysm3 = imp2ss(B_IRF(:,10)*tSIRF,tSIRF);

% Test Result
[y_New3,t_New3]=impulse(sysm3);

%%
figure; 
plot(B_IRFt(:,thetaIndx),B_IRF(:,thetaIndx)); hold on;
plot(t_New3,y_New3)