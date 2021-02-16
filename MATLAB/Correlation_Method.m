% ME 155C Control System Lab Project: Non-Parametric Identification Correlation Method
% By: Alex Nguyen

clc; clear; close all;

load('Part1.mat')

%FREQUENCY DATA
item = sort([10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 2 30 35 ...
    3 40 45 4 50 55 5 60 65 6 70 75 7 80 8 90 9]);

%% Correlation Method
alpha = 2; %input amplitude
H = zeros(length(item),2); %preallocation

for i = 1:length(item)
    %DEFINE INPUT & OUTPUT DATA
    u = eval(sprintf('Part1_%drad_s_input_1',item(i))); %input [V]
    y = encoder(eval(sprintf('Part1_%drad_s_output_1',item(i)))); %output [s,m,m]
    t = y(:,1); %time vector [s]
    w = item(i); %continuous time frequency [rad/s]
    
    %     t = y(1:end-1,1); %time vector [s]
    %     y = y(2:end,2:3) - y(1:end-1,2:3); %assuming integrator at z = 1
    
    %ADJUST DATA
    [~,ia] = findpeaks(u); %local max
    n1 = ia(1); %new index starting point
    n2 = ia(end); %index ending point
    
    %CORRECTING PHASE IN DATA
    u = u(n1:n2); y = y(n1:n2,2:3); t = t(n1:n2);
    
%     u = u(n1:n2); y = y(n1:n2,:); t = t(n1:n2); %assuming integrator at z = 1 
    
    %CROP THE DATA
    T = t(end)-t(1); %input sine-wave duration [s]
    Ts = t(2)-t(1); %sample time [s]
    Omega = Ts*w; %discrete time angular velocity [rad]
    f = w/2/pi; %frequency [Hz] - specifically 10 rad/s
    tpeak = floor((T-3)*f)/f; %keep last 3 s of experiment
    kpeak = round(tpeak/Ts); %convert to samples
        
    %CORRELATION METHOD
    ycrop = y(kpeak:end,1:2); %cropped output data [m,m]
    k = 1:length(ycrop); %sample number
    z = exp(-1i*Omega*k);
    K = 1/length(ycrop)*(z*ycrop);
    H(i,:) = 2*K/alpha; %H estimate
end

%CONSTRUCTING BODE PLOT
Hdb = 20*log10(abs(H)); %Gain [dB]
Hdeg1 = phase(H(:,1))*180/pi; %Cart x1 Phase [deg]
Hdeg2 = phase(H(:,2))*180/pi; %Cart x2 Phase [deg]
w = sort(item); %frequency values [rad/s]

%IDEAL TRANSFER FUNCTION
G = tf(2.97*61.2,[1 13.24 127.15 810.37 0]);
[mag,ang,wout] = bode(G); %bode mag & phase values
mag = squeeze(mag); phase = squeeze(ang);
mag = 20*log10(abs(mag)); 

figure;
subplot(2,1,1)
semilogx(w,Hdb(:,2),wout,mag,'--r'); grid on;
xlabel('w [rad/s]')
ylabel('Magnitude [dB]')
legend('Estimated (x_2)','Ideal TF','location','best')
xlim([2 100])
subplot(2,1,2)
semilogx(w,Hdeg2,wout,phase,'--r'); grid on;
xlabel('w [rad/s]')
ylabel('Phase [deg]')
legend('Estimated (x_2)','Ideal TF', 'location','best')
xlim([2 100])
sgtitle('Bode Plot - Spring Cart (x_2) Comparison') 

%% ESTIMATED TRANSFER FUNCTION 
gain = 10.^(Hdb(:,2)/20);
response = gain.*exp(1i*Hdeg2*pi/180);
gfr = idfrd(response,w,Ts);
x2_est = tfest(gfr,4,0); %estimated x2 c-t transfer function

%ESTIMATE TRANSFER FUNCTION ZPK
b = [366.5]; %numerator coefficients
a = [1 11.65 487.8 3216 121.7]; %denominator coefficients
[z1,p1,k1] = tf2zp(b,a);
sys_est1 = zpk(z1,p1,k1);
save('Process1.mat','sys_est1','z1','p1','k1','Ts') %saving .mat data

%ESTIMATED BODE MAGNITUDE AND PHASE DATA
[m,a,wnew] = bode(x2_est); %estimated x2 bode values
m = squeeze(m); a = squeeze(a);
m = 20*log10(abs(m)); 

%BODE PLOT COMPARISON
figure; 
subplot(2,1,1)
semilogx(wnew,m,wout,mag,'--r'); grid on;
xlabel('w [rad/s]')
ylabel('Magnitude [dB]')
legend('Approximate TF','Ideal','location','best')
xlim([2 100])
subplot(2,1,2)
semilogx(wnew,a,wout,phase,'--r'); grid on;
xlabel('w [rad/s]')
ylabel('Phase [deg]')
legend('Approximate TF','Ideal', 'location','best')
xlim([2 100])
sgtitle('Bode Plot - Estimated TF (x_2) Comparison') 