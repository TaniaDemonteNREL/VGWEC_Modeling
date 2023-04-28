%% Model Reduction of Radiation State Space

thetaN = 10; % Discretized Theta Increments in degrees
thetaVector = 0:dTheta:90; 
ssOrder     = 6; 
tfRadAR_RO  = cell(1,10); 

for i = 1:length(thetaVector)
 
    sys1 = ss(tfRadAR{i},tfRadBR{i},tfRadCR{i},tfRadDR{i});  

    % Reduce the order of the state space model using balanced truncation
    sys_reduced = balred(sys1,ssOrder);

    tfRadAR_RO{i} = sys_reduced.A;  
    tfRadBR_RO{i} = sys_reduced.B; 
    tfRadCR_RO{i} = sys_reduced.C; 
    tfRadDR_RO{i} = sys_reduced.D; 

    figure; bode(sys1,'r',sys_reduced,'b');
    legend('Original','Reduced');

end 

%% Save Data
save('flap2_Frad_Reduced','tfRadAR_RO','tfRadBR_RO','tfRadCR_RO','tfRadDR_RO'); 

%% Format Data
tfRadAR_ROm = cell2mat(tfRadAR_RO); 
tfRadBR_ROm = cell2mat(tfRadBR_RO); 
tfRadCR_ROm = cell2mat(tfRadCR_RO); 
tfRadDR_ROm = cell2mat(tfRadDR_RO); 

% Compare the frequency responses of the original and reduced models
% figure; bode(sys1,'r',sys_reduced,'b');
% legend('Original','Reduced');