close all; clear all; clc;

M = 0.5;
m = 0.2;
b = 0.1;
I = 0.006;
g = 9.8;
L = 0.3;

q= (M+m)*(I+m*L^2)-(m*L)^2;
s = tf('s');
P_pend= (m*L*s/q)/(s^3 + (b*(I + m*L^2))*s^2/q - ((M + m)*m*g*L)*s/q - b*m*g*L/q);
Kp = 100;
Ki = 1;
Kd = 70;

C = pid(Kp, Ki, Kd);
T = feedback(P_pend, C);

t = 0:0.01:10;

figure(1);
impulse(T,t);
axis([0, 2.5, -0.2, 0.2]);
grid on;
title({'Response of Pendulum Position to an Impulse Disturbance'; ...
    'under PID Control: Kp = 100, Ki = 1, Kd = 70'})

figure(2)
rlocus(T)