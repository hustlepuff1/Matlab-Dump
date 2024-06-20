close all; clear all; clc;

M = 0.5;    %mass of the cart
m = 0.2;    %mass of the pendulum
b = 0.1;    %friction of cart
I = 0.006;  %inertia of pendulum
g = 9.8;
L = 0.3;
q= (M+m)*(I+m*L^2)-(m*L)^2;
s = tf('s');

P_cart = (((I+m*L^2)/q)*s^2 - (m*g*L/q))/(s^4 + (b*(I + m*L^2))*s^3/q - ((M + m)*m*g*L)*s^2/q - b*m*g*L*s/q);

P_pend= (m*L*s/q)/(s^3 + (b*(I + m*L^2))*s^2/q - ((M + m)*m*g*L)*s/q - b*m*g*L/q);

sys_tf = [P_cart; P_pend];

inputs = {'u'};
outputs = {'x'; 'phi'};

set(sys_tf, 'InputName', inputs)
set(sys_tf,'OutputName', outputs)

sys_tf

t = 0:0.01:1;
figure;
impulse(sys_tf, t);
grid on;
title('Open-Loop Impulse Response')

[zeros poles] = zpkdata(P_pend, 'v')
[zeros1 poles1] = zpkdata(P_cart, 'v')

