%%
clc; close all;

m = 1;
k = 1;
c = 0.2;
F = 1;

A = [0 1; -k/m -c/m];
B = [0 1/m]';
C = [1 0];
D = [0];

sys = ss(A,B,C,D);

s = tf('s');
sys = 1/(m*s^2+c*s+k);

%%
clc; close all;

s = tf('s');
P = 1/(s^2+10*s+20);
step(P)

%%
clc; close all;

s = tf('s');
P = 1/(s^2+10*s+20);

Kp = 300;
C = pid(Kp);
T = feedback(C*P,1);

t = 0:0.01:2;
step(T,t)

%%
clc; close all;

s = tf('s');
P = 1/(s^2+10*s+20);

Kp = 300;
Kd = 10;
C = pid(Kp,0,Kd);
T = feedback(C*P,1);

t = 0:0.01:2;
step(T,t)

%%
clc; close all;

s = tf('s');
P = 1/(s^2+10*s+20);

Kp = 30;
Ki = 70;
C = pid(Kp,Ki);
T = feedback(C*P,1);

t = 0:0.01:2;
step(T,t)

%%
clc; close all;
s = tf('s');
H = (s+2)/(s^2 + 2*s +3);
z_ol = zero(H);
p_ol = pole(H);
axis([-10 2 -2 2])
K = 0.1;
T = feedback(K*H, 1);
p_cl = pole(T);
realpart = real(p_cl);
imagpart = imag(p_cl);
grid on
hold on
plot(real(z_ol), imag(z_ol), 'bO')
plot(real(p_ol), imag(p_ol), 'bX')

figure(2)
rlocus(H)
axis([-10 2 -2 2])

%%
clc; close all;

s = tf('s');
sys = (s + 7)/(s*(s + 5)*(s + 15)*(s + 20));
rlocus(sys)
axis([-22 3 -15 15])
zeta = 0.7;
wn = 1.8;
sgrid(zeta,wn)

[k, poles] = rlocfind(sys);

K = 350;
sys_cl = feedback(K*sys,1);
step(sys_cl)

%%
clc; close all;
c = 2;
k = 8;
m = 2;