%% Interpolation of Time-Varying Hydrodynamic Coefficients
% This script will interpolate the hydrodynamic coefficients for each
% discrete angle. 

%% Theta Interpolation

currentTheta = 5:dTheta:95; 
fineTheta    = currentTheta(1):new_dTheta:currentTheta(end); 
LW = 'Linewidth';
%% Added Mass at Infinite Frequency

AinfInterp = interp1(currentTheta,A_PitchInf,fineTheta);
[AinfFit, S, a_mu] = polyfit(currentTheta, A_PitchInf, polyOrder); 

hydroVCoeff.AinfInterp = AinfInterp;
hydroVCoeff.AinfFit    = AinfFit; 
hydroVCoeff.a_mu       = a_mu;

figure; subplot(2,1,1); plot(fineTheta,polyval(AinfFit, fineTheta,[], a_mu), LW, 1.5);
hold on;  plot(fineTheta,AinfInterp,'-.',currentTheta,A_PitchInf,'-o', LW, 1.5)
title(' Added Mass Moment of Inertia at Inf Freq'); grid on;
xlabel('Flap Angle (deg)'); ylabel('Added mass (kg-m^2/rad)');
legend('PolyFit','Interpolation','Actual Discrete Data')

%% Hydrostatic Coefficient 
c       = hydroVCoeff.c; 
cInterp = interp1(currentTheta,c,fineTheta); 
[cFit, ~, c_mu] = polyfit(currentTheta,c,polyOrder); 

hydroVCoeff.cInterp = cInterp;
hydroVCoeff.cFit    = cFit; 
hydroVCoeff.c_mu    = c_mu;


subplot(2,1,2);
plot(fineTheta,polyval(cFit, fineTheta,[], c_mu), LW, 1.5);
hold on;  plot(fineTheta,cInterp,'-.',currentTheta,c,'-o', LW, 1.5)
title(' Hydrostatic Coefficient'); grid on;
xlabel('Flap Angle (deg)'); ylabel('Hydrostatic coefficient (N-m/rad)');
legend('PolyFit','Interpolation','Actual Discrete Data')


%% Radiation Damping Coefficient
 
bRadInterp = interp1(currentTheta,B1_all',fineTheta);
hydroVCoeff.bRadInterp = bRadInterp'; 


figure; subplot(2,1,1)
h1 = plot(w,bRadInterp'./1000,'-.',LW, 1.5); hold on; 
h2 = plot(w,B1_all./1000,'-*', LW, 2);
title('Pitch Radiation Damping Coefficient'); grid on;
xlabel('Flap Angle (degrees)'); ylabel('Radiation Damping (kN-m/(rad/s))');
legend([h1(1) h2(1)],'Interp1','Actual Data')
xlim([0,15])
%% Excitation Moment

TExc       = hydroVCoeff.KE; 
FExcInterp = interp1(currentTheta,TExc',fineTheta);
hydroVCoeff.FExcInterp = FExcInterp'; 

subplot(2,1,2)
plot(w,FExcInterp'./1000,'-.',LW, 1.5); hold on; 
plot(w,TExc./1000,'-*', LW, 2)
title('Pitch Excitation Moment Magnitude'); grid on;
xlabel('Flap Angle (degrees)'); ylabel('Excitation Moment (kN-m)');
legend([h1(1) h2(1)],'Interp1','Actual Data')
xlim([0,15])

%% Radiation Damping 

% First, we get the transfer function coefficients from the SS models. 
for i = 1:length(currentTheta)    
    [bR(i,:), aR(i,:)] = ss2tf(tfRadAR_RO{i},tfRadBR_RO{i},tfRadCR_RO{i},tfRadDR_RO{i});
end 

%% Then, we interpolate between the coefficients. 

for i = 1:size(aR,2)

    aRInterp(i,:) = interp1(currentTheta,aR(:,i),fineTheta);
    [aRFit(i,:), ~, aR_mu(i,:)] = polyfit(currentTheta, aR(:,i), polyOrder); 

    bRInterp(i,:) = interp1(currentTheta,bR(:,i),fineTheta);
    [bRFit(i,:), ~, bR_mu(i,:)] = polyfit(currentTheta, bR(:,i), polyOrder); 

    figure; subplot(2,1,1); plot(fineTheta,polyval(aRFit(i,:), fineTheta,[], aR_mu(i,:)), LW, 1.5);
    hold on;  plot(fineTheta,aRInterp(i,:),'-.',currentTheta,aR(:,i),'o', LW, 1.5)
    title("aR " + num2str(i) + " Coefficient"); grid on;
    legend('PolyFit','Interp1','Actual Data'); xlabel('Pitch Angle (degrees)')

    subplot(2,1,2);  plot(fineTheta,polyval(bRFit(i,:), fineTheta,[], bR_mu(i,:)), LW, 1.5);
    hold on;  plot(fineTheta,bRInterp(i,:),'-.',currentTheta,bR(:,i),'o', LW, 1.5)
    title("bR " + num2str(i) + " Coefficient"); grid on;
    legend('PolyFit','Interp1','Actual Data'); xlabel('Pitch Angle (degrees)'); 
    
end

%%
hydroVCoeff.aRFit = aRFit; 
hydroVCoeff.bRFit = bRFit; 
hydroVCoeff.aR_mu = aR_mu; 
hydroVCoeff.bR_mu = bR_mu; 
hydroCoeff.fineTheta  = fineTheta; 