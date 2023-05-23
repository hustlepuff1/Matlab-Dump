clc; close all;

num = 8;
den = conv([1 2], [1 8]);
oltf = tf(num, den);

rlocus(oltf); grid

k = 3.44;
cltf = feedback(k*oltf,1);

figure; step(cltf); grid