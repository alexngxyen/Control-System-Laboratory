% ME 155C Control System Lab Project: Parametric Identification Chirp &
% Square Wave
% By: Alex Nguyen

clc; clear; close all;

load('Part2.mat')

%% Ideal Transfer Function
G = tf(2.97*61.2,[1 13.24 127.15 810.37 0]);
[m,ang] = bode(G); %bode mag & phase values
m = 20*log10(abs(squeeze(m))); phase = squeeze(ang);

%% SQUARE WAVE - INPUT/OUTPUT DATA
u1 = Part2_square_dif_10_input_1; %input
y1 = encoder(Part2_square_dif_10_output_1); %output
u2 = Part2_square_10_input_1; %input
y2 = encoder(Part2_square_10_output_1); %output

%TIME
t = y2(:,1); %time vector [s]
Ts = t(2)-t(1); %sample time [s]

% ADJUST DATA - DEALING WITH KNOWN PARAMETERS
%assume the process has an integrator (d-t pole at z = 1) in its TF
y1 = y1(:,3); y2 = y2(:,3); %redefining data
ybar1 = y1(2:end) - y1(1:end-1); ybar2 = y2(2:end) - y2(1:end-1); %output [m] 
ubar1 = u1(1:end-1); ubar2 = u2(1:end-1); %input [V]

% ADJUST DATA - SCALING OUTPUT SIGNAL
alpy = 7500;
alpu = 1000;
ybar1 = ybar1*alpy; ybar2 = ybar2*alpy; %scaled output
ubar1 = ubar1*alpu; ubar2 = ubar2*alpu; %scaled input

% ADJUST DATA - DOWN SAMPLING
L = round(0.04/Ts); %down sample from 1kHz to 40 Hz 
ybar1 = resample(ybar1,1,L,1); ybar2 = resample(ybar2,1,L,1); %downsampling output [V]
ubar1 = resample(ubar1,1,L,1); ubar2 = resample(ubar2,1,L,1); %downsampling input [V]
Tnew = L*Ts;

%% DEVELOP ARX MODEL FOR DIFFERENT SIGNALS OF X2 (CART ATTATCHED TO SPRING)
%CREATE THE BEST ARX MODEL 
na = 3; nb = 4; nk = 0; 

dat = {0}; %preallocation
dat{1} = iddata(ybar1,ubar1,Tnew);
dat{2} = iddata(ybar2,ubar2,Tnew);

data = merge(dat{1},dat{2});
model = arx(data,[na nb nk]); %estimated model
z = tf('z',Tnew); %used to create a d-t TF model 
sysd = alpy/alpu/(z-1)*tf(model,'Measured'); %d-t transfer function
P0_square = d2c(sysd,'zoh'); %c-t transfer function

%% CHIRP SIGNAL -INPUT/OUTPUT DATA
% INPUT
u = zeros(9001,6); %preallocation

%STORING INPUT DATA
u(:,1) = Part2_chirp_0p001_1_rad_s_input_1; %chirp signal [0.001 Hz,1 Hz]
u(:,2) = Part2_chirp_0p001_25_rad_s_input_1; %chirp signal [0.001 Hz,25 Hz]
u(:,3) = Part2_chirp_0p001_90_input_1; %chirp signal [0.001 Hz,90 Hz]
u(:,4) = Part2_chirp_10_50_input_1; %chirp signal [10 Hz,50 Hz]
u(:,5) = Part2_chirp_1_10_rad_s_input_1; %chirp signal [1 Hz,10 Hz]
u(:,6) =  Part2_chirp_20_90_input_1; %chirp signal [20 Hz,90 Hz]

%OUTPUT
y = zeros(9001,12); %preallocation

%STORING OUTPUT DATA - UNITS [s,m,m]
y(:,1) = Part2_chirp_0p001_1_rad_s_output_1(:,3); %chirp signal [0.001 Hz,1 Hz]
y(:,2) = Part2_chirp_0p001_25_rad_s_output_1(:,3); %chirp signal [0.001 Hz,25 Hz]
y(:,3) = Part2_chirp_0p001_90_output_1(:,3); %chirp signal [0.001 Hz,90 Hz]
y(:,4) = Part2_chirp_10_50_output_1(:,3); %chirp signal [10 Hz,50 Hz]
y(:,5) = Part2_chirp_1_10_rad_s_output_1(:,3); %chirp signal [1 Hz,10 Hz]
y(:,6) =  Part2_chirp_20_90_output_1(:,3); %chirp signal [20 Hz,90 Hz]

%CONVERT TICKS TO METERS
encoder = 0.1/4096; %[rev/ticks]
y = encoder*y; 

% ADJUST DATA - DEALING WITH KNOWN PARAMETERS
%assume the process has an integrator (d-t pole at z = 1) in its TF
ybar = zeros(length(y)-1,size(y,2)); ubar = ybar;
for i = 1:6
    ybar(:,i) = y(2:end,i) - y(1:end-1,i); %output
    ubar(:,i) = u(1:end-1,i); %input
end

% ADJUST DATA - SCALING OUTPUT SIGNAL
alpy = 7500;
alpu = 1000;
ybar = alpy.*ybar(:,1:6); %scaled output
ubar = alpu.*ubar(:,1:6); %scaled input

% ADJUST DATA - DOWN SAMPLING
L = round(0.04/Ts); %down sample from 1kHz to 40 Hz 
Tnew = L*Ts; %new sample time [s]
y = zeros(length(resample(ybar(:,1),1,L,1)),6); u = y; %preallocation
for i = 1:size(ybar,2)
    y(:,i) = resample(ybar(:,i),1,L,1); %output resampled
    u(:,i) = resample(ubar(:,1),1,L,1); %input resampled
end

%% DEVELOP ARX MODEL FOR DIFFERENT SIGNALS OF X2 (CART ATTATCHED TO SPRING)
%CREATE THE BEST ARX MODEL 
na = 3; nb = 4; nk = 0; 

dat = {0}; %preallocation
for i = 1:6
    dat{i} = iddata(y(:,i),u(:,i),Tnew); 
end
data = merge(dat{1},dat{2},dat{3},dat{4},dat{5},dat{6});

model = arx(data,[na nb nk]); %estimated model
z = tf('z',Tnew); %used to create a d-t TF model 
sysd = alpy/alpu/(z-1)*tf(model,'Measured'); %d-t transfer function
P0_chirp = d2c(sysd,'zoh'); %c-t transfer function

%PLOT BODE 
w = logspace(-1,4,100); %frequency [rad/s]
[m1,a1] = bode(G,w); m1 = 20*log10(abs(squeeze(m1))); a1 = squeeze(a1); %ideal
[m2,a2] = bode(P0_chirp,w); m2 = 20*log10(abs(squeeze(m2))); a2 = squeeze(a2) - 360; %parametric - chirp
[m3,a3] = bode(P0_square,w); m3 = 20*log10(abs(squeeze(m3))); a3 = squeeze(a3) - 360; %parametric - square

figure;
subplot(2,1,1)
semilogx(w,m1,w,m2,w,m3); grid on;
ylabel('Magnitude [dB]')
xlabel('frequency [rad/s]')
subplot(2,1,2)
semilogx(w,a1,w,a2,w,a3); grid on;
ylabel('Phase [deg]')
xlabel('frequency [rad/s]')
legend('Ideal','Parametric - chirp','Parametric - square','Location','best')
sgtitle('Bode Plot - Parametric Identification')

%SAVE .MAT DATA
save('ParametricTF.mat','P0_chirp','P0_square','G') %saving .mat data