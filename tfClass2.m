clc; close all;
num = 2;
den = 1;
sysG1 = tf(num, den);

num = 4;
den = [1 0];
sysG2 = tf(num, den);

sysG3 = parallel(sysG1, sysG2);

num = 1;
dem = [1 0];
sysG4 = tf(num, den);

sysG5 = series(sysG3, sysG4);

num = 1;
den = 1;
sysG6 = tf(num, den);

sysCL = feedback(sysG5, sysG6, -1);
