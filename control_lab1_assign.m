clc; close all;

m = 2;
c = 2;
k = 8;
s = tf('s');
sys = 1/(m*s^2 + c*s + k);

% num = 1;
% den = [2 2 8];
% sys = tf(num, den);

Kp = 100;
Ki = 80;
Kd = 35;
K = pid(Kp, Ki, Kd);
T = feedback(K*sys,1);

t = 0:0.1:8;

figure(1)
step(T,t)
stepinfo(T)

figure(2)
rlocus(T)