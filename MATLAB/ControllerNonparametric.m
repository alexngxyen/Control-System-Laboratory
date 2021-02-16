% ME 155C Control System Lab Project: Controller Design
% By: Alex Nguyen

clc; clear; close all;

%NONPARAMETRIC AND PARAMETRIC PROCESS TRANSFER FUNCTION
load('Process1.mat'); P0a = sys_est1; %nonparametric process
load('Process2.mat'); P0b = sys2; %parametric process - all data
load('ParametricTF.mat'); %parametric process - chirp and square data TF

%IDEAL TRANSFER FUNCTION
s = tf('s'); %ct variable 's'
G = tf(2.97*61.2,[1 13.24 127.15 810.37 0]); 

%REQUIRED PHASE MARGIN 
OS = 15; %percent overshoot
zeta = -log(OS/100)/sqrt(pi^2+log(OS/100)^2); %damping ratio
pm_req = atand(2*zeta/sqrt(-2*zeta^2+sqrt(1+4*zeta^4))); %phase margin required

%DEVLOPING CONTROLLER - LEAD COMPENSATOR
K = 24; %proportional gain 
[~,Pm,~,~] = margin(K*P0a); %phase margin
phi = pm_req - Pm + 15; %maximum phase margin
alpha = (1-sind(phi))/(1+sind(phi));
[~,~,~,wm] = margin(K*P0a/sqrt(alpha)); %new crossover frequency [rad/s]
Kc = K/alpha;
T = 1/wm/sqrt(alpha); 
C = Kc*(alpha*T*s+1)/(s*T+1); %lead compensator

%STEP RESPONSE INFO
L = C*P0a; %open loop gain
v = .85; %voltage step input [V]
Y1 = v*feedback(L,1); %closed loop - controller 
Y2 = v*feedback(G,1); %closed loop feedback - no controller
stepinfo(Y1) %step-response characterisitics

figure;
step(Y1,Y2); %step-response output - Estimated & Ideal
legend('Controller','No Controller','location','best')
title('Closed-Loop Step Response with Ideal TF')

%ROBUSTNESS WITH RESPECT TO MEASUREMENT NOISE
N1 = v*feedback(C*P0a,1); %nonparametric
N2 = v*feedback(C*P0_square,1); %parametric - square signal
Y = v*feedback(C*G,1); %ideal 

b = [1416 3618];
a = [0.3819 6.056 61.79 436.6 2476 4256];
[z,p,k] = tf2zp(b,a);
sys = zpk(z,p,k);

figure;
step(N1,N2,Y2)
legend('Nonpar','Par - Square','Ideal','location','best')

%SAVE .MAT DATA
save('controller.mat','C')