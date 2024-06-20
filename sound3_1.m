clc; clear all; close all;

t = [0:0.001:1];
x = [];
x_tmp = zeros(size(t));

for n = 1:2:39
    x_tmp = x_tmp + 4/pi*(1/n*sin(2*pi*n*t));
    x = [x;
        x_tmp];
end

plot(t, x(3,:), t, x(7,:), t,x(20,:))
xlabel('\itt\rm(seconds)');
ylabel('\itx\rm(\itt\rm)')