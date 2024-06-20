clear all; clc; close all;

Ac=1; Am=4; fm=1; fc=8;
t=0:0.0001:4;
x=Ac*cos(2*pi*fc*t + Am*sin(2*pi*fm*t));
a=hilbert(x);
fx=diff(unwrap(angle(a)))./diff(t)/(2*pi);

subplot(3,1,1)
plot(t, abs(a), t, x, 'r:'); axis([0 4 -1.5 1.5])
xlabel('Time (s)'); ylabel('\itA x\rm(\itt\rm)')

subplot(3,1,2)
plot(t, unwrap(angle(a))); axis([0 4 0 220])
xlabel('Time (s)');
ylabel('\it\phi x\rm(\itt\rm)')

subplot(3,1,3)
plot(t(2:end),fx); axis([0 4 0 13])
xlabel('Time (s)'); ylabel('\itf x\rm(\itt\rm)')
