%% 1번 예제
clc; close all;
figure(1);
num = 1;
den = conv([1 1 0], [1 5]);
oltf = tf(num, den);

rlocus(oltf); grid
%% 2번 예제
clc; close all;
figure(1);
num = conv([1 1], [1 0 81]);
den = conv([1 13], [1 0 100 0 0]);
oltf = tf(num, den);

rlocus(oltf); grid

figure(2);
k1 = 89.1;
k2 = 32.4;
cltf1 = feedback(k1*oltf,1);
cltf2 = feedback(k2*oltf,1);
hold on; grid on;
step(cltf1)
step(cltf2)
legend('k = 89.1', 'k = 32.4');
%% 3번 예제
clc; close all;
figure(1);
num = 1;
den = conv([1 -1], [1 2]);
oltf = tf(num, den);

rlocus(oltf); grid

figure(2);
k = 2.5;
cltf = feedback(k*oltf,1);
step(cltf);