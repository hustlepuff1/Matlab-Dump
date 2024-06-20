close all; clear all; clc;

M = 0.5;    %mass of the cart
m = 0.2;    %mass of the pendulum
b = 0.1;    %friction of cart
I = 0.006;  %inertia of pendulum
g = 9.8;
L = 0.3;

q= (M+m)*(I+m*L^2)-(m*L)^2;
s = tf('s');
P_pend= (m*L*s/q)/(s^3 + (b*(I + m*L^2))*s^2/q - ((M + m)*m*g*L)*s/q - b*m*g*L/q);
Kp = 100;
Ki = 1;
Kd = 20;

C = pid(Kp, Ki, Kd);
T = feedback(P_pend, C);

t = 0:0.01:10;

figure(1);
impulse(T,t);
axis([0, 2.5, -0.2, 0.2]);
grid on;
title({'Response of Pendulum Position to an Impulse Disturbance'; ...
    'under PID Control: Kp = 100, Ki = 1, Kd = 20'})

figure(2);
rlocus(P_pend)
grid on;
title('Root Locus of Plant (under Proportional Control')

figure(3);
z = [-3 -4];
p = 0;
k = 1;
C = zpk(z, p, k);
rlocus(C*P_pend)
grid on;
title('Root Locus with PID Controller')

figure(4);
[k, poles] = rlocus(P_pend);
K = 20;
T = feedback(P_pend, K*C);
impulse(T)
title('Impulse Disturbance Response of Pendulum Angle under PID Control');

% figure(5);
% P_cart = (((I+m*L^2)/q)*s^2 - (m*g*L/q))/(s^4 + (b*(I + m*L^2))*s^3/q - ((M + m)*m*g*L)*s^2/q - b*m*g*L*s/q);
% T2 = feedback(1, P_pend*C)*P_cart;
% t2 = 0:0.01:8.5;
% impulse(T2, t2);
% title('Impulse Disturbance Response of Cart Position under PID Control');