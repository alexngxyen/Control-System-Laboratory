% ME 155C Control System Lab Project: Closed-Loop Performace
% By: Alex Nguyen

clc; clear; close all;

%LOADING DATA
load('controller.mat') %lead compensator controller
load('Process1.mat') %nonparametric identification
load('Process2.mat') %parametric identification - all data
load('ParametricTF.mat') %parametric identification - chirp and square data

%EVALUATE CLOSED LOOP STEP RESPONSE - ADJUSTED (v = 0.85)
v = 0.85; %input voltage
Y1 = v*feedback(G*C,1); %ideal closed loop
Y2 = v*feedback(sys_est1*C,1); %nonparametric closed loop
Y3 = v*feedback(P0_square*C,1); %parametric closed loop - square input

ideal = stepinfo(Y1); %ideal
nonparametric_step_info = stepinfo(Y2); %nonparametric
parametric_step_info = stepinfo(Y3); %parametric - square input

figure;
step(Y3);
title(sprintf('Closed-Loop Step-Response with Input %4.2f V',v))

%CLOSED-LOOP FREQUENCY RESPONSE
w = logspace(-1,3,100); %frequency [rad/s]
[m1,a1] = bode(Y1,w); m1 = 20*log10(abs(squeeze(m1))); a1 = squeeze(a1); %ideal
[m2,a2] = bode(Y2,w); m2 = 20*log10(abs(squeeze(m2))); a2 = squeeze(a2); %nonparametric
[m3,a3] = bode(Y3,w); m3 = 20*log10(abs(squeeze(m3))); a3 = squeeze(a3) - 360; %parametric

figure;
subplot(2,1,1)
semilogx(w,m3); grid on;
ylabel('Magnitude [dB]')
xlabel('frequency [rad/s]')
subplot(2,1,2)
semilogx(w,a3); grid on;
ylabel('Phase [deg]')
xlabel('frequency [rad/s]')
sgtitle('CL Frequency Response - Parametric with Controller')

figure;
subplot(2,1,1)
semilogx(w,m1,w,m2,w,m3); grid on;
ylabel('Magnitude [dB]')
xlabel('frequency [rad/s]')
subplot(2,1,2)
semilogx(w,a1,w,a2,w,a3); grid on;
ylabel('Phase [deg]')
xlabel('frequency [rad/s]')
legend('Ideal','Nonparametric','Parametric','Location','best')
sgtitle('Closed-Loop Frequency Response')

%PRINT RESULTS
nonparametric_step_info 
parametric_step_info