clear; clc; close all;
x = 0:0.1:10;
y = exp(-0.5.*x).*(0.16*cos(0.5*sqrt(35).*x)+0.027*sin(0.5*sqrt(35).*x));
plot(x, real(y));