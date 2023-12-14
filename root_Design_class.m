clc; close all;

z1 = 1;
p1 = 12;
num1 = [1 z1];
den1 = [1 p1 0 0];
oltf1 = tf(num1, den1);

z2 = 1;
p2 = 4;
num2 = [1 z2];
den2 = [1 p2 0 0];
oltf2 = tf(num2, den2);

z3 = 1;
p3 = 9;
num3 = [1 z3];
den3 = [1 p3 0 0];
oltf3 = tf(num3, den3);

figure(1);
rlocus(oltf1); hold on
rlocus(oltf2);
rlocus(oltf3);
axis equal
grid on

figure(2);
k1 = 78;
cltf1 = feedback(k1*oltf1, 1);
k2 = 24.7;
cltf2 = feedback(k2*oltf2, 1);
k3 = 40;
cltf3 = feedback(k3*oltf3, 1);
step(cltf1); hold on
step(cltf2);
step(cltf3);





