% ME 155C Control System Lab Project: Controller Design
% By: Alex Nguyen

clc; clear; close all;

%NONPARAMETRIC AND PARAMETRIC PROCESS TRANSFER FUNCTION
load('Process1.mat'); P0a = zpk(z1,p1,k1); %nonparametric process
load('Process2.mat'); P0b = zpk(z2,p2,k2); %parametric process - all data
load('ParametricTF.mat'); %parametric process - chirp and square data TF

%IDEAL TRANSFER FUNCTION
s = tf('s'); %ct variable 's'
G = tf(2.97*61.2,[1 13.24 127.15 810.37 0]); 

%REQUIRED PHASE MARGIN 
OS = 15; %percent overshoot
zeta = -log(OS/100)/sqrt(pi^2+log(OS/100)^2); %damping ratio
pm_req = atand(2*zeta/sqrt(-2*zeta^2+sqrt(1+4*zeta^4))); %phase margin required

%DEVLOPING CONTROLLER - LEAD COMPENSATOR
K = 20; %proportional gain 
[~,Pm,~,~] = margin(K*G); %phase margin
phi = pm_req - Pm + 10; %maximum phase margin
alpha = (1-sind(phi))/(1+sind(phi));
[~,~,~,wm] = margin(K*G/sqrt(alpha)); %new crossover frequency [rad/s]
T = 1/wm/sqrt(alpha); 
C = 12*(alpha*T*s+1)/(s*T+1); %lead compensator 

%STEP RESPONSE INFO
L = C*G; %open loop gain
v = .85; %voltage step input [V]
Y1 = v*feedback(L,1); %closed loop - controller 
Y2 = v*feedback(G,1); %closed loop feedback - no controller
stepinfo(Y1) %step-response characterisitics

figure;
step(Y1,Y2); %step-response output - Estimated & Ideal
legend('Controller','No Controller','location','best')
title('Closed-Loop Step Response with Ideal TF')

% %GAIN & PHASE MARGIN
% [Gm,Pm,Wcg,Wcp] = margin(L); 
% fprintf('Gain Margin: %4.4f, Phase Margin: %4.4f\n',Gm,Pm)
% fprintf('Gain Crossover Freq: %4.4f, Phase Crossover Freq: %4.4f\n',Wcg,Wcp)

%ROBUSTNESS WITH RESPECT TO MEASUREMENT NOISE
Y2 = v*feedback(C*P0_square,1); %parametric - square signal
Y3 = v*feedback(C*P0a,1); %nonparametric

figure;
step(Y3,Y2,Y1)
legend('Nonpar','Par - Square','Ideal','location','best')

%SAVE .MAT FILE
save('controller.mat','C')