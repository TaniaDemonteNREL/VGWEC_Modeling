function [w, E]  = getExcitationMagnitude(varargin,plotFlag,rho,g)
% Plots the excitation force magnitude for each hydro structure's bodies in
% the heave, surge and pitch degrees of freedom.
% 
% Usage:
% ``plotExcitationMagnitude(hydro, hydro2, hydro3, ...)``
% 
% Parameters
% ----------
%     varargin : struct(s)
%         The hydroData structure(s) created by the other BEMIO functions.
%         One or more may be input.
% 

if isempty(varargin)
    error(['plotExcitationMagnitude: No arguments passed. Include one or more hydro ' ...
        'structures when calling: plotExcitationMagnitude(hydro1, hydro2, ...)']);
end

B=1;  % Wave heading index

if plotFlag
figHandle = figure('Position',[950,500,975,521]);
titleString = ['Normalized Excitation Force Magnitude: ',...
    '$$\bar{X_i}(\omega,\theta) = {\frac{X_i(\omega,\theta)}{{\rho}g}}$$'];
subtitleString = {'Surge','Heave','Pitch'};
xString = {'$$\omega (rad/s)$$','$$\omega (rad/s)$$','$$\omega (rad/s)$$'};
yString = {['$$\bar{X_1}(\omega,\theta$$',' = ',num2str(varargin{1}.theta(B)),'$$^{\circ}$$)'],...
    ['$$\bar{X_3}(\omega,\theta$$',' = ',num2str(varargin{1}.theta(B)),'$$^{\circ}$$)'],...
    ['$$\bar{X_5}(\omega,\theta$$',' = ',num2str(varargin{1}.theta(B)),'$$^{\circ}$$)']};

notes = {''};
end 

numHydro = length(varargin);
for ii = 1:numHydro
    numBod = varargin(ii).Nb;
    tmp1 = strcat('X',num2str(ii));
    X.(tmp1) = varargin(ii).w;
    tmp2 = strcat('Y',num2str(ii));
    a = 0;
    for i = 1:numBod
        m = varargin(ii).dof(i);
        Y.(tmp2)(1,i,:) = squeeze(varargin(ii).ex_ma(a+1,B,:));
        Y.(tmp2)(2,i,:) = squeeze(varargin(ii).ex_ma(a+3,B,:));
        Y.(tmp2)(3,i,:) = squeeze(varargin(ii).ex_ma(a+5,B,:));
        legendStrings{i,ii} = [varargin(ii).body{i}];
        a = a + m;
    end
end

if plotFlag
formatPlot(figHandle,titleString,subtitleString,xString,yString,X,Y,legendStrings,notes);
saveas(figHandle,'Excitation_Magnitude.png');
end 

w = X.X1; 

E.Surge = squeeze(Y.Y1(1,1,:)) * rho * g;
E.Heave = squeeze(Y.Y1(2,1,:)) * rho * g;
E.Pitch = squeeze(Y.Y1(3,1,:)) * rho * g;

end
