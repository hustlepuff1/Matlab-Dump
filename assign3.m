clc; close all;

%% 3.1
k = 2;

num = k;
den = [1 2 k];
sysA = tf(num, den);

t = 0:0.01:10;

step(sysA, t)
stepinfo(sysA)

%% 3.2
k = 2.44;

num = k;
den = [1 2 k];
sysA = tf(num, den);

t = 0:0.01:10;

step(sysA, t)
stepinfo(sysA)
