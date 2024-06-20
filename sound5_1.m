clear all; clc; close all;

Ac=1; Am=0.5; fm=1; fc=10;
t=0:0.001:3;
x=(Ac+Am*cos(2*pi*fm*t)).*cos(2*pi*fc*t);
a=hilbert(x);
fx=diff(unwrap(angle(a)))./diff(t)/(2*pi);

subplot(3,1,1)
plot(t, abs(a), t, x, 'r:')
axis([0 3 -2 2])
xlabel('Time (s)'); ylabel('\itA x\rm(\itt\rm)')

subplot(3,1,2)
plot(t, unwrap(angle(a)))
axis([0 3 0 200])
xlabel('Time (s)'); ylabel('\it\phi x\rm(\itt\rm)')

subplot(3,1,3)
plot(t(2:end),fx,'r')
axis([0 3 8 12])
xlabel('Time (s)'); ylabel('\itf x\rm(\itt\rm)')