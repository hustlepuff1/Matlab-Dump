clc; close all;
num = 2;
den = [1 3 2];
sys = tf(num, den);
[z, p, k] = tf2zp(num, den);
t = 0:0.01:10;
impulse(sys, t);