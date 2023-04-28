%% Radiation Damping 
clear bR2 aR2   aRFit aR_mu bRFit bR_mu

bR2 = zeros(10,7); 
aR2 = zeros(10,7); 

% First, we get the transfer function coefficients from the SS models. 
for i = 1:length(currentTheta)    
    n = length(radA{i}); 

    [bR2(i,end-(n):end), aR2(i,end-(n):end)] = ss2tf(radA{i},radB{i},radC{i},0);
end 

%% Then, we interpolate between the coefficients. 

polyOrder = 10; 

for i = 1:size(aR2,2)

    aRInterp2(i,:) = interp1(currentTheta,aR2(:,i),fineTheta);
%     [aRFit(i,:), ~, aR_mu(i,:)] = polyfit(currentTheta, aR(:,i), polyOrder); 

    bRInterp2(i,:) = interp1(currentTheta,bR2(:,i),fineTheta);
%     [bRFit(i,:), ~, bR_mu(i,:)] = polyfit(currentTheta, bR(:,i), polyOrder); 

    figure; subplot(2,1,1); %plot(fineTheta,polyval(aRFit(i,:), fineTheta,[], aR_mu(i,:)), LW, 1.5);
    hold on;  plot(fineTheta,aRInterp2(i,:),'-.',currentTheta,aR2(:,i),'o', LW, 1.5)
    title("aR " + num2str(i) + " Coefficient"); grid on;
    legend('Interp1','Actual Data')

    subplot(2,1,2);  %plot(fineTheta,polyval(bRFit(i,:), fineTheta,[], bR_mu(i,:)), LW, 1.5);
    hold on;  plot(fineTheta,bRInterp2(i,:),'-.',currentTheta,bR2(:,i),'o', LW, 1.5)
    title("bR " + num2str(i) + " Coefficient"); grid on;
    legend('Interp1','Actual Data')

    
end

%%
hydroVCoeff.aRFit = aRFit; 
hydroVCoeff.bRFit = bRFit; 
hydroVCoeff.aR_mu = aR_mu; 
hydroVCoeff.bR_mu = bR_mu;

%% New rad SS from SS fitting toolbox
load('C:\Users\tdemonte\OneDrive - NREL\VGWEC\Full Flap Data\Flap2_FDI_Frad.mat')

AR = radA{2};  
BR = radB{2}; 
CR = radC{2}; 
DR = zeros(1,1);

n = 4 + 1;
AR2 = zeros(6,6); 
AR2(end-n+2:end,end-n+2:end) = AR; 

BR2 = zeros(6,1); 
BR2(end-n+2:end,1) =  BR; 

CR2 = zeros(1,6);
CR2(1,end-n+2:end) = CR; 

DR2 = DR;

%%
radCoeff = hydroVCoeff.b(w_WaveIndx,2)*1000; 