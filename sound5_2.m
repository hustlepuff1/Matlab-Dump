clear all; clc; close all;

% parameter
A=1;
zeta=0.01;
fn=10;
wn=2*pi*fn;
wd=wn*sqrt(1-zeta^2);
phi=0;

t=0:0.001:6;

x = A*exp(-zeta*wn*t).*sin(wd*t+phi);

a=hilbert(x);
ax=log(abs(a));

subplot(2,1,1)
plot(t, abs(a), t, x, 'r:'); axis([0 6 -1.5 1.5])
xlabel('Time (s)');
ylabel('\itA x\rm(\itt\rm)')

subplot(2,1,2)
plot(t, ax); axis([0 6 -6 1])
xlabel('Time (s)');
ylabel('ln\itA x\rm(\itt\rm)')

p=polyfit(t(1000:4000), ax(1000:4000), 1);

format long
zeta_est=-p(1)/wn