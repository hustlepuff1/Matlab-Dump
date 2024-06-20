clc; clear all; close all;

L_del_a = 2.0; L_p = -0.5;

% inner loop
num = L_del_a;
den = [1 -L_p];
oltf1 = tf(num, den);

k_as = 6.82;

cltf1 = feedback(k_as * oltf1, 1);
oltf2 = cltf1 * tf(1, [1 0]);

k = 7.33;
cltf2 = feedback(k * oltf2, 1);


figure(1);
rlocus(oltf1)
figure(2);
rlocus(oltf2)
figure(3);
rlocus(cltf2)
figure(4);
step(cltf2)


