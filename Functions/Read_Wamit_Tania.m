%% Read Wamit and Frequency Domain
% Read the coefficients from WAMIT Outputs
% Files needed: fname.3; fname.1; fname.hst
% Last modified by:  Sal Husain
% Last updated on: November 1st, 2017


% clc
% clear

%%
tic

filename = 'flap1_0half';
fname = filename;

%% EXTRACT VOLUME
fileID = fopen([fname '.mmx']);

mmx_file = textscan(fileID,'%s','HeaderLines',1);
GRAV     = str2double(mmx_file{1, 1}{2, 1});
ULEN     = str2double(mmx_file{1, 1}{5, 1});
VOLX     = str2double(mmx_file{1, 1}{24, 1});
VOLY     = str2double(mmx_file{1, 1}{25, 1});
VOLZ     = str2double(mmx_file{1, 1}{26, 1});

RHO      = 1000; % Changed from 1000
RAD = 2.5;
fclose(fileID);
%% Read excitation force - Diffraction Potential 
disp('Importing excitation coefficients')

n_angle = 1;
n_body = 1;    %number of bodies
modes_tmp = [1 2 3 4 5 6];
modes = [];
for i = 1:n_body
    modes = [modes modes_tmp];
end
n_modes = length(modes);

% % n_freq = 161;
mod_vec_denorm_full = [2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3];
mod_vec_denorm = mod_vec_denorm_full(modes);
fileID = fopen([fname '.3']);
wam_exc_out = textscan(fileID,'%f %f %f %f %f %f %f','HeaderLines',1);
fclose(fileID);
wam_exc.data = cell2mat(wam_exc_out);

n_modes = find( (wam_exc.data(2:end,3)==wam_exc.data(1,3)),1 );

I = wam_exc.data(1:n_modes,3);
T = wam_exc.data(1:n_modes:end,1);
Ww = 2.*pi./T;

n_freq = size(T,1);
E_tmp = zeros(6,1,n_freq);
Phi_tmp = E_tmp;
REpart_temp = E_tmp;
IMpart_tmp = E_tmp;

for ii = 1:n_modes
    ind = ii ;
    E_tmp(I(ii),:) = wam_exc.data(ind:n_modes:end,4);
    Phi_tmp(I(ii),:) = wam_exc.data(ind:n_modes:end,5);
    REpart_temp(I(ii),:) = wam_exc.data(ind:n_modes:end,6);
    IMpart_tmp(I(ii),:) = wam_exc.data(ind:n_modes:end,7);
end


E = RHO * GRAV * E_tmp;
Phi = Phi_tmp;
REpart = RHO * GRAV * REpart_temp;
IMpart = RHO * GRAV * IMpart_tmp;
toc

%% Read Froude-Krylov component of the Excitation Force - Diffraction Potential

disp('Importing non linear excitation force  coefficients')
% read excitation
% wam_exc = importdata([fname '.3']);
n_angle = 1;
n_body = 1;    %number of bodies
modes_tmp = [1 2 3 4 5 6];
modes = [];
for i = 1:n_body
    modes = [modes modes_tmp];
end
n_modes = length(modes);

% % n_freq = 161;
mod_vec_denorm_full = [2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3];
mod_vec_denorm = mod_vec_denorm_full(modes);
fileID = fopen([fname '.3fk']);
wam_exc_out = textscan(fileID,'%f %f %f %f %f %f %f','HeaderLines',1);
fclose(fileID);
wam_exc.data = cell2mat(wam_exc_out);

n_modes = find( (wam_exc.data(2:end,3)==wam_exc.data(1,3)),1 );

I = wam_exc.data(1:n_modes,3);
T = wam_exc.data(1:n_modes:end,1);
Ww = 2*pi/T;

n_freq = size(T,1);
E_tmpfk = zeros(6,1,n_freq);
Phi_tmpfk = E_tmpfk;
REpart_tempfk = E_tmpfk;
IMpart_tmpfk = E_tmpfk;

for ii = 1:n_modes
    ind = ii ;
    E_tmpfk(I(ii),:) = wam_exc.data(ind:n_modes:end,4);
    Phi_tmpfk(I(ii),:) = wam_exc.data(ind:n_modes:end,5);
    REpart_tempfk(I(ii),:) = wam_exc.data(ind:n_modes:end,6);
    IMpart_tmpfk(I(ii),:) = wam_exc.data(ind:n_modes:end,7);
end


E_fk = RHO * GRAV * E_tmpfk;
Phi_fk = Phi_tmpfk;
REpart_fk = RHO * GRAV * REpart_tempfk;
IMpart_fk = RHO * GRAV * IMpart_tmpfk;
E1_fk = squeeze(E_fk(3,1,:));
toc



%% Read Scattering component of the Excitation Force - Diffraction Potential

disp('Importing non linear excitation force  coefficients')
% read excitation
% wam_exc = importdata([fname '.3']);
n_angle = 1;
n_body = 1;    %number of bodies
modes_tmp = [1 2 3 4 5 6];
modes = [];
for i = 1:n_body
    modes = [modes modes_tmp];
end
n_modes = length(modes);

% % n_freq = 161;
mod_vec_denorm_full = [2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3];
mod_vec_denorm = mod_vec_denorm_full(modes);
fileID = fopen([fname '.3sc']);
wam_exc_out = textscan(fileID,'%f %f %f %f %f %f %f','HeaderLines',1);
fclose(fileID);
wam_exc.data = cell2mat(wam_exc_out);

n_modes = find( (wam_exc.data(2:end,3)==wam_exc.data(1,3)),1 );

I = wam_exc.data(1:n_modes,3);
T = wam_exc.data(1:n_modes:end,1);
Ww = 2*pi/T;

n_freq = size(T,1);
E_tmpsc = zeros(6,1,n_freq);
Phi_tmpsc = E_tmpsc;
REpart_tempsc = E_tmpsc;
IMpart_tmpsc = E_tmpsc;

for ii = 1:n_modes
    ind = ii ;
    E_tmpsc(I(ii),:) = wam_exc.data(ind:n_modes:end,4);
    Phi_tmpsc(I(ii),:) = wam_exc.data(ind:n_modes:end,5);
    REpart_tempsc(I(ii),:) = wam_exc.data(ind:n_modes:end,6);
    IMpart_tmpsc(I(ii),:) = wam_exc.data(ind:n_modes:end,7);
end


E_sc = RHO * GRAV * E_tmpsc;
Phi_sc = Phi_tmpsc;
REpart_sc = RHO * GRAV * REpart_tempsc;
IMpart_sc = RHO * GRAV * IMpart_tmpsc;
E1_sc = squeeze(E_sc(3,1,:));
toc

%% Read Scattering component of the Excitation Force - Haskind 
% 
% disp('Importing non linear excitation force  coefficients')
% % read excitation
% % wam_exc = importdata([fname '.3']);
% n_angle = 1;
% n_body = 1;    %number of bodies
% modes_tmp = [1 2 3 4 5 6];
% modes = [];
% for i = 1:n_body
%     modes = [modes modes_tmp];
% end
% n_modes = length(modes);
% 
% % % n_freq = 161;
% mod_vec_denorm_full = [2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3];
% mod_vec_denorm = mod_vec_denorm_full(modes);
% fileID = fopen([fname '.2sc']);
% wam_exc_out = textscan(fileID,'%f %f %f %f %f %f %f','HeaderLines',1);
% fclose(fileID);
% wam_exc.data = cell2mat(wam_exc_out);
% 
% n_modes = find( (wam_exc.data(2:end,3)==wam_exc.data(1,3)),1 );
% 
% I = wam_exc.data(1:n_modes,3);
% T = wam_exc.data(1:n_modes:end,1);
% Ww = 2*pi/T;
% 
% n_freq = size(T,1);
% E_tmpsc = zeros(6,1,n_freq);
% Phi_tmpsc = E_tmpsc;
% REpart_tempsc = E_tmpsc;
% IMpart_tmpsc = E_tmpsc;
% 
% for ii = 1:n_modes
%     ind = ii ;
%     E_tmpsc(I(ii),:) = wam_exc.data(ind:n_modes:end,4);
%     Phi_tmpsc(I(ii),:) = wam_exc.data(ind:n_modes:end,5);
%     REpart_tempsc(I(ii),:) = wam_exc.data(ind:n_modes:end,6);
%     IMpart_tmpsc(I(ii),:) = wam_exc.data(ind:n_modes:end,7);
% end
% 
% 
% E_scH = RHO * GRAV * E_tmpsc;
% Phi_sc = Phi_tmpsc;
% REpart_sc = RHO * GRAV * REpart_tempsc;
% IMpart_sc = RHO * GRAV * IMpart_tmpsc;
% E1_scH = squeeze(E_scH(3,1,:));
% toc

%% Read Froude-Krylov component of the Excitation Force - Haskind
% 
% disp('Importing non linear excitation force  coefficients')
% % read excitation
% % wam_exc = importdata([fname '.3']);
% n_angle = 1;
% n_body = 1;    %number of bodies
% modes_tmp = [1 2 3 4 5 6];
% modes = [];
% for i = 1:n_body
%     modes = [modes modes_tmp];
% end
% n_modes = length(modes);
% 
% % % n_freq = 161;
% mod_vec_denorm_full = [2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3];
% mod_vec_denorm = mod_vec_denorm_full(modes);
% fileID = fopen([fname '.2fk']);
% wam_exc_out = textscan(fileID,'%f %f %f %f %f %f %f','HeaderLines',1);
% fclose(fileID);
% wam_exc.data = cell2mat(wam_exc_out);
% 
% n_modes = find( (wam_exc.data(2:end,3)==wam_exc.data(1,3)),1 );
% 
% I = wam_exc.data(1:n_modes,3);
% T = wam_exc.data(1:n_modes:end,1);
% Ww = 2.*pi./T;
% 
% n_freq = size(T,1);
% E_tmpfk = zeros(6,1,n_freq);
% Phi_tmpfk = E_tmpfk;
% REpart_tempfk = E_tmpfk;
% IMpart_tmpfk = E_tmpfk;
% 
% for ii = 1:n_modes
%     ind = ii ;
%     E_tmpfk(I(ii),:) = wam_exc.data(ind:n_modes:end,4);
%     Phi_tmpfk(I(ii),:) = wam_exc.data(ind:n_modes:end,5);
%     REpart_tempfk(I(ii),:) = wam_exc.data(ind:n_modes:end,6);
%     IMpart_tmpfk(I(ii),:) = wam_exc.data(ind:n_modes:end,7);
% end
% 
% 
% E_fkH = RHO * GRAV * E_tmpfk;
% Phi_fk = Phi_tmpfk;
% REpart_fk = RHO * GRAV * REpart_tempfk;
% IMpart_fk = RHO * GRAV * IMpart_tmpfk;
% E1_fkH = squeeze(E_fkH(3,1,:));
% toc

%% Read excitation force - Haskind
% disp('Importing excitation coefficients')
% % read excitation
% % wam_exc = importdata([fname '.3']);
% n_angle = 1;
% n_body = 1;    %number of bodies
% modes_tmp = [1 2 3 4 5 6];
% modes = [];
% for i = 1:n_body
%     modes = [modes modes_tmp];
% end
% n_modes = length(modes);
% 
% % % n_freq = 161;
% mod_vec_denorm_full = [2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3];
% mod_vec_denorm = mod_vec_denorm_full(modes);
% fileID = fopen([fname '.2']);
% wam_exc_out = textscan(fileID,'%f %f %f %f %f %f %f','HeaderLines',1);
% fclose(fileID);
% wam_exc.data = cell2mat(wam_exc_out);
% 
% n_modes = find( (wam_exc.data(2:end,3)==wam_exc.data(1,3)),1 );
% 
% I = wam_exc.data(1:n_modes,3);
% T = wam_exc.data(1:n_modes:end,1);
% Ww = 2*pi/T;
% 
% n_freq = size(T,1);
% E_tmp = zeros(6,1,n_freq);
% Phi_tmp = E_tmp;
% REpart_temp = E_tmp;
% IMpart_tmp = E_tmp;
% 
% for ii = 1:n_modes
%     ind = ii ;
%     E_tmp(I(ii),:) = wam_exc.data(ind:n_modes:end,4);
%     Phi_tmp(I(ii),:) = wam_exc.data(ind:n_modes:end,5);
%     REpart_temp(I(ii),:) = wam_exc.data(ind:n_modes:end,6);
%     IMpart_tmp(I(ii),:) = wam_exc.data(ind:n_modes:end,7);
% end
% 
% 
% EH = RHO * GRAV * E_tmp;
% Phi = Phi_tmp;
% REpart = RHO * GRAV * REpart_temp;
% IMpart = RHO * GRAV * IMpart_tmp;
% E1H = squeeze(EH(3,1,:));
% toc

%% Read Response Amplitude Operator - Haskind
% disp('Importing Response Amplitude Operator coefficients')
% 
% n_angle = 1;
% n_body = 1;    %number of bodies
% modes_tmp = [1 2 3 4 5 6];
% modes = [];
% for i = 1:n_body
%     modes = [modes modes_tmp];
% end
% n_modes = length(modes);
% 
% % % n_freq = 161;
% mod_vec_denorm_full = [2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3 2 2 2 3 3 3];
% mod_vec_denorm = mod_vec_denorm_full(modes);
% fileID = fopen([fname '.4']);
% wam_exc_out = textscan(fileID,'%f %f %f %f %f %f %f','HeaderLines',1);
% fclose(fileID);
% wam_exc.data = cell2mat(wam_exc_out);
% 
% n_modes = find( (wam_exc.data(2:end,3)==wam_exc.data(1,3)),1 );
% 
% I = wam_exc.data(1:n_modes,3);
% T = wam_exc.data(1:n_modes:end,1);
% Ww = 2*pi/T;
% 
% n_freq = size(T,1);
% RAO_tmp = zeros(6,1,n_freq);
% Phi_tmp = RAO_tmp;
% REpart_temp = RAO_tmp;
% IMpart_tmp = RAO_tmp;
% 
% for ii = 1:n_modes
%     ind = ii ;
%     RAO_tmp(I(ii),:) = wam_exc.data(ind:n_modes:end,4);
%     Phi_tmp(I(ii),:) = wam_exc.data(ind:n_modes:end,5);
%     REpart_temp(I(ii),:) = wam_exc.data(ind:n_modes:end,6);
%     IMpart_tmp(I(ii),:) = wam_exc.data(ind:n_modes:end,7);
% end
% 
% 
% RAO = RAO_tmp; %*RHO * GRAV * 
% Phi = Phi_tmp;
% REpart =  REpart_temp; %*RHO * GRAV * 
% IMpart =  IMpart_tmp; %*RHO * GRAV * 
% RAO1 = squeeze(RAO(3,1,:));
% toc
%% Read Hydrodynamic Pressure 
% disp('Importing Hydrodynamic Pressure Coefficients')
% 
% 
% n_body = 1;    %number of bodies
% modes_tmp = [1 2 3 4 5 6];
% modes = [];
% for i = 1:n_body
%     modes = [modes modes_tmp];
% end
% n_modes = length(modes);
% 
% 
% fileID = fopen([fname '.5p']);
% wam_exc_out = textscan(fileID,'%f %f %f %f %f %f %f %f','HeaderLines',1);
% fclose(fileID);
% wam_exc.data = cell2mat(wam_exc_out);
% 
% n_modes = 1521; % number of panels
% 
% I = wam_exc.data(1:n_modes,3);
% T = wam_exc.data(1:n_modes:end,1);
% Ww = 2*pi/T;
% 
% n_freq = size(T,1);
% E_tmp = zeros(6,1,n_freq);
% Phi_tmp = E_tmp;
% REpart_temp = E_tmp;
% IMpart_tmp = E_tmp;
% 
% for ii = 1:n_modes
%     ind = ii ;
%     E_tmp(I(ii),:) = wam_exc.data(ind:n_modes:end,5);
%     Phi_tmp(I(ii),:) = wam_exc.data(ind:n_modes:end,6);
%     REpart_temp(I(ii),:) = wam_exc.data(ind:n_modes:end,7);
%     IMpart_tmp(I(ii),:) = wam_exc.data(ind:n_modes:end,8);
% end
% 
% 
% P.Mod = RHO * GRAV * E_tmp;
% P.Phi = Phi_tmp;
% P.REpart = RHO * GRAV * REpart_temp;
% P.IMpart = RHO * GRAV * IMpart_tmp;
% toc
%% read radiation

disp('Importing radiation coefficients ')


fileID = fopen([fname '.1']);
wam_rad_out = textscan(fileID,'%f %f %f %f %f','HeaderLines',1);
fclose(fileID);

wam_rad.data = cell2mat(wam_rad_out);

n_modes = find(and(wam_rad.data(2:end,2)==wam_rad.data(1,2), wam_rad.data(2:end,3)==wam_rad.data(1,3)),1);

I = wam_rad.data(1:n_modes,2);
J = wam_rad.data(1:n_modes,3);

I0 = find(wam_rad.data(:,1)>0,1,'first');

T = wam_rad.data(I0:n_modes:end,1);
w = 2*pi./T;

n_freq = size(wam_rad.data(I0:end,:),1)/n_modes;
A_tmp = zeros(6, 6, n_freq);
B_tmp = A_tmp;

% multiply damping coefficients by OMEGA - for denormalization
%B_vect = 2*pi * wam_rad.data(:,5) ./ wam_rad.data(:,1);

for ii = 1:n_modes
    ind = ii + I0-1;
    A_tmp(I(ii),J(ii),:) = wam_rad.data(ind:n_modes:end,4);
    B_tmp(I(ii),J(ii),:) = wam_rad.data(ind:n_modes:end,5);    
end

norm_mat = RHO*ULEN.^kron([3,4;4,5], ones(3)); 

% de-normalization matrix: RHO * L^k 

%A_tmp = A_tmp*RHO;
%B_tmp = B_tmp*RHO;
toc

%% read added mass for 0 and Inf frequency


A0_tmp = zeros(6,6);
Ainf_tmp = A0_tmp;


for ii = 1:I0-1
    i = wam_rad.data(ii,2);
    j = wam_rad.data(ii,3);
    if wam_rad.data(ii,1) < 0
        A0_tmp(i,j) = wam_rad.data(ii,4);
    else
        Ainf_tmp(i,j) = wam_rad.data(ii,4);
    end
end

A0 = A0_tmp *RHO;
Ainf = Ainf_tmp *RHO;
A0 = A0_tmp;
Ainf = Ainf_tmp;
toc
%% read hydrostatic restoring coefficients from hst
disp('Hydrostatic restoring coefficients ')

fileID = fopen([fname '.hst']);
Hydrostatic_out = textscan(fileID,'%f %f %f','HeaderLines',1);
fclose(fileID);
wam_hydro.data = cell2mat(Hydrostatic_out);
n_modes = (6*n_body)^2;
I = wam_hydro.data(1:n_modes,1);
J = wam_hydro.data(1:n_modes,2);
for ii = 1:n_modes
    ind = ii;
    Kh_tmp(I(ii),J(ii),:) = wam_hydro.data(ind:n_modes:end,3);    
end

KH = RHO * GRAV * Kh_tmp;                          
%denormalize

toc
%%

freq = w;
freq_exc = w;
% A = cat(3, zeros(6*n_body,6*n_body), A_tmp);
% B = cat(3, zeros(6*n_body,6*n_body), B_tmp);
% save('Testsal.mat', 'A', 'Ainf', 'B', 'E','Phi', 'freq','freq_exc','KH')

A = A_tmp;
B = B_tmp;
toc

disp(' ')
disp('done')
disp(' ')
disp('*********************** ')


%%
close all
figure(1)
A1 = squeeze(A(3,3,:));
% A2 = squeeze(A(9,9,:));
% A3 = squeeze(A(15,15,:));
%A1 = A1 * (2/3) * pi * (1)^3; 
hold on
plot(freq,A1,'+r')%,freq,A2,'+b',freq,A3,'ok')/(RHO*pi*8*2/3)
title('Added mass')
ylabel('Added mass (kg)')
xlabel('\omega rad/s')
legend('body 1','body 2','body 3')
grid on
saveas(gcf,strcat(filename,'_A1','.fig'))


%%
figure(2)
B1 = squeeze(B(3,3,:));
% B2 = squeeze(B(9,9,:));
% B3 = squeeze(B(15,15,:));
%B1 = B1 * (2/3) *pi * (1)^3;
plot(freq,B1,'ok')%,freq,B2,'ob',freq,B3,'-k')/(RHO*pi*8*2/3)
title('Radiation Damping')
ylabel('Radiation Damping (Ns/m)')
xlabel('\omega rad/s')
legend('body 1','body 2','body 3')
grid on
saveas(gcf,strcat(filename,'_B1','.fig'))

%%
figure(3)
E1 = squeeze(E(3,1,:));
% E2 = squeeze(E(9,1,:));
% E3 = squeeze(E(15,1,:));

plot(freq,E1,'+b')%,freq_exc,E2,'ob',freq_exc,E3,'-k')/(RHO*pi*8*2/3)
title('Excitation Force Coefficient(N/m)')
xlabel('\omega rad/s')
ylabel('Exciting Force Coefficient (N/m)')
legend('body 1','body 2','body 3')
grid on
saveas(gcf,strcat(filename,'_E','.fig'))

%%
figure(4)
plot(freq,E1_fk,'+b',freq,E1_sc,'or')%,freq_exc,E2,'ob',freq_exc,E3,'-k')/(RHO*pi*8*2/3)
title('FK and Scattering Components of Excitation Force (N/m)')
xlabel('\omega rad/s')
ylabel('Exciting Force Coefficient (N/m)')
legend('Fexc Froude Krylov Component','Fexc Scattering Component')
grid on
saveas(gcf,strcat(filename,'_Efk_AND_Esc','.fig'))

%%

figure(5)
plot((RAD/9.81)*freq.^2,A1/((VOLX+VOLY+VOLZ)*RHO),'+r')
hold on
Bs = abs(B1./(abs(freq)*(VOLX+VOLY+VOLZ)*RHO));
plot((RAD/9.81)*freq.^2,Bs,'ok')
legend('Added Mass coefficient k ','Radiation Damping coefficient e')
xlabel ('\omega^2 a/g')
grid on

%%

figure(6)

plot(freq,A1,'+r')

hold on

plot(freq,B1,'ok')

legend('Added Mass Kg ','Radiation Damping N/(m/s)')
xlabel ('\omega rads/s')
grid on

saveas(gcf,strcat(filename,'_HD_coeffs','.fig'))
%%
% figure(7)
% P1 = squeeze(P.Mod(1,1,:));
% P1 = unique(P1,'stable');
% % B2 = squeeze(B(9,9,:));
% % B3 = squeeze(B(15,15,:));
% %B1 = B1 * (2/3) *pi * (1)^3;
% plot(freq,P1)%,freq,B2,'ob',freq,B3,'-k')/(RHO*pi*8*2/3)
% title('Hydrodynamic Pressure')
% ylabel('Hydrodynamic Pressure')
% xlabel('\omega rad/s')
% %legend('body 1','body 2','body 3')
% grid on
% saveas(gcf,strcat(filename,'_P1','.fig'))

%%
try
save(filename,'A','B','A1','B1','E','E1',...
   'KH','Phi','REpart','T','w');
catch
  save(filename,'A','B','E','E_fk','E_sc'...
   ,'KH','Phi','Phi_fk','Phi_sc','REpart','REpart_fk','REpart_sc','T','w');  
end