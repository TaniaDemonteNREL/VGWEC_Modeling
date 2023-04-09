%% Convolution of Rad and Exc Forces. 
load('Flap2_DataALL.mat');

t        = 0:0.01:50;
eta      = amp_Wave * sin(w_Wave*t); 
velocity = 10 * cos(w_Wave*t); 

thetaIndx = 10; 
B_irf     = real(B_IRF(:,thetaIndx));
E_irf     = real(E_IRF(end/2+1:end,thetaIndx));


%% Convolution 1
y1R       = zeros(1,length(B_irf)+length(velocity)-1);
y1E       = zeros(1,length(E_irf)+length(eta)-1);

for i = 1 : length(velocity)
for j = 1 : length(B_irf)    
    
y1R(i+j-1)  = y1R(i+j-1) +  B_irf(j)*velocity(i); 

end
end

Frad_1 = y1R(1,:);

for i = 1 : length(eta)
for j = 1 : length(E_irf)     
y1E(i+j-1)  = y1E(i+j-1) +  E_irf(j)*eta(i); 
end
end

Fexc_1 = y1E(1,:);
%% Convolution 2l
y2R = conv(B_irf,velocity);
y2E = conv(E_irf,eta);

%% Convolution 3
y3R = conv(B_Pitch(w_WaveIndx,thetaIndx),velocity);
y3E = conv(E_Pitch(w_WaveIndx,thetaIndx),eta);

%% Convolution 4
y4R = B_Pitch(w_WaveIndx,thetaIndx) * velocity;
y4E = E_Pitch(w_WaveIndx,thetaIndx) * eta;

%% Figure RADIATION

LW = 'Linewidth'; 

figure; subplot(2,1,1); plot(t,y1R(1:length(t))./1000,LW,1.5); hold on; 
plot(t,y2R(1:length(t))./1000,'-.',LW,1.5); legend('HandMade Conv','MATLAB Conv')
ylabel('Rad Force (kNm)'); xlabel('Time (s)'); grid on;
sgtitle('Radiation Force')

subplot(2,1,2);plot(t,y3R(1:length(t)),LW,1.5); hold on;
plot(t,y4R(1:length(t)),'-.',LW,1.5);  grid on;
legend('MATLAB Conv - Single Freq','Multiplication')
ylabel('Rad Force (Nm)'); xlabel('Time (s)'); 

%% Figure EXCITATION

LW = 'Linewidth'; 

figure; subplot(2,1,1); plot(t,y1E(1:length(t))./1000,LW,1.5); hold on; 
plot(t,y2E(1:length(t))./1000,'-.',LW,1.5); legend('HandMade Conv','MATLAB Conv')
ylabel('Exc Force (kNm)'); xlabel('Time (s)'); grid on;
sgtitle('Excitation Force')

subplot(2,1,2);plot(t,y3E(1:length(t)),LW,1.5); hold on;
plot(t,y4E(1:length(t)),'-.',LW,1.5);  grid on;
legend('MATLAB Conv - Single Freq','Multiplication')
ylabel('Rad Force (Nm)'); xlabel('Time (s)'); 