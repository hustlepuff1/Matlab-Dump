clear all; clc; close all;
fs = 500;
t=-5:1/fs:5;

lambda = 300;
t0 = 0.15;
a = 0.2;

x=exp(-lambda*abs(t));
y = x + a*exp(-lambda*abs(t-t0));

X=fft(x);
Y=fft(y);

N = length(x);
fp = 0:fs/N:fs/2;
fn = -fs/N:-fs/N:-fs/2;
f = [fliplr(fn) fp];

plot(f, fftshift(abs(X)/fs), 'r.')
xlabel("Frequency (Hz)")
ylabel("Modulus")

hold on

plot(f, fftshift(abs(Y)/fs))
hold off


