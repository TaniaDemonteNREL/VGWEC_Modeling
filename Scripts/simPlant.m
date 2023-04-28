function thetaD = simPlant(x0,t,dt,w_Wave, hydroCoeff, hydroVCoeff, Fexc, Frad, pFlag)
  
  % number of points to simulate
  N = length(t);
  
  % grab memory for the results
  xSim = zeros(N,length(x0));
  
  % install the initial condition
  xSim(1,:) = x0;
  
  f = @vgWecPlant; % this gives just a tad of flexibility...
  
  % simulate using ode3
  for i = 1:length(t)-1
    xn   = xSim(i,:);    
    fn   = f( xn, w_Wave, hydroCoeff, hydroVCoeff, Fexc(i), Frad(i));
    fnh  = f( xn+dt*fn/2   , w_Wave, hydroCoeff, hydroVCoeff, Fexc(i), Frad(i));
    fnhh = f( xn+3*dt*fnh/4, w_Wave, hydroCoeff, hydroVCoeff, Fexc(i), Frad(i));
  
    xSim(i+1,:) = xn + dt/9*(2*fn + 3*fnh + 4*fnhh);
  end
  
 % Make a plot if requested. This should not be used during optimization
   % as it would be very slow.
   if pFlag == 1
     plot(t,xSim.*180/pi,'LineWidth',1.5); grid on; 
     legend('theta','thetaDot'); xlabel('Time'); ylabel('Degrees / Degrees/s');

   end
   
   %thetaD = simPlant(theta00,t,dt,w_Wave, hydroCoeff, hydroVCoeff, FE, FR, 1);
   thetaD = xSim; 
end 