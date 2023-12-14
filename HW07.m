%% 문제 1번
clc; close all;

num = [100 100];
den = [1 110 1000];
oltf = tf(num, den);
margin(oltf);
grid on

%% 문제 2번
clc; close all;

num = [10 100];
den = conv([1 -1 0], [0.01 1]);
oltf = tf(num, den);
figure(1);
margin(oltf);
grid on

figure(2);
rlocus(oltf);
grid on


