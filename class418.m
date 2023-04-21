clc; clear all; close all;
numG = 8;
denG = [1 10 16];
G = tf(numG, denG);
t = 0:0.1:5;

Kp = 1:2:9;
for i = 1:length(Kp)
    numD = Kp(i);
    denD = 1;
    D = tf(numD, denD);
    cltf = feedback(D*G, 1);
    y(:, i) = step(cltf, t);
end

plot(t, y, 'LineWidth', 2);


