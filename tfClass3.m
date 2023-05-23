clc; close all;
wn = 2;
zeta = 0:0.1:1;
t = 0:0.1:5;

num = wn^2;
for i = 1:length(zeta);
    den = [1 2*zeta(i)*wn wn^2];
    sys = tf(num, den);
    step(sys, t);
    hold on; pause

 
end