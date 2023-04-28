LW = 'LineWidth'; 
constantData = amp_2_K_5_w_2; 
T_V_Data     = amp_2_K_5_w_2_TV;
%% Figure 1

figure; 
subplot(2,1,1); plot(constantData{12}.Values,'b',LW,1.5); 
hold on; plot(T_V_Data{12}.Values,'r',LW,1.5); 
xlim([200,300]); title('Total Force'); grid on; 
legend('Single Hyrdrodynamics','T-V Hydrodynamics')
xlabel('Time (s)'); ylabel('Force (N)');

subplot(2,1,2); plot(constantData{11}.Values,'b',LW,1.5); 
hold on; plot(T_V_Data{11}.Values,'r',LW,1.5); 
xlim([200,300]); title('Flap Rotation'); grid on; 
xlabel('Time (s)'); ylabel('Angular Displacement (deg)')

%% Figure 2

figure; 
subplot(2,1,1); plot(constantData{3}.Values,'b',LW,1.5); 
hold on; plot(T_V_Data{3}.Values,'r',LW,1.5); 
xlim([200,300]); title('Excitation Force'); grid on; 
xlabel('Time (s)'); ylabel('Force (N)');
legend('Single Hyrdrodynamics','T-V Hydrodynamics'); 

subplot(2,1,2); plot(T_V_Data{11}.Values,'r',LW,1.5); 
xlim([200,300]); title('Flap Rotation'); grid on; 
xlabel('Time (s)'); ylabel('Angular Displacement (deg)')
%% Figure 3

figure; subplot(2,1,1); plot(constantData{4}.Values,'b',LW,1.5); 
hold on; plot(T_V_Data{4}.Values,'r',LW,1.5); 
xlim([200,300]); title('Radiation Force'); grid on; 
xlabel('Time (s)'); ylabel('Force (N)')

subplot(2,1,2); plot(T_V_Data{11}.Values,'r',LW,1.5); 
xlim([200,300]); title('Flap Rotation'); grid on; 
xlabel('Time (s)'); ylabel('Angular Displacement (deg)')