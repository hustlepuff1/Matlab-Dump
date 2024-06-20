clc; close all; clear;

s = tf('s');
P_pitch = (1.151*s+0.1774) / (s^3 + 0.739*s^2+0.921*s);
sys_cl = feedback(P_pitch,1);
t = 0:0.01:10;

figure(1)
step(0.2*P_pitch,t);
axis([0 10 0 0.8]);
ylabel('Pitch angle (rad)');
title("Open Loop Step Response");

A = [-0.313 56.7 0; 
-0.0139 -0.426 0; 
0 56.7 0];
B = [0.232; 0.0203; 0];
C = [0 0 1];
D = 0;
pitch_ss = ss(A,B,C,D);

figure(2)
step(0.2*pitch_ss,t)
axis([0 10 0 0.8]);
ylabel('pitch angle (rad)')
title("Openloop Step Response")

figure(3);
step(0.2*sys_cl);