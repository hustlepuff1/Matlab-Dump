clc; clear all; close all;

fs = 500;
T = 10;
t = 0:1/fs:T-1/fs;
p = 40;

x = sin(2*pi*p*t);

x1 = x(1:5:end);
x2 = x(1:10:end);

N = length(x);
N1 = length(x1);
N2 = length(x2);

X = fft(x);
X1 = fft(x1);
X2 = fft(x2);

f = fs*(0:N-1)/N;
f1 = (fs/5)*(0:N1-1)/N1;
f2 = (fs/10)*(0:N2-1)/N2;

subplot(3,1,1);
plot(f, abs(X)/N)
xlabel('Frequency (Hz)');
ylabel('Modulus')
axis([0 500 0 0.55])
subplot(3,1,2);
plot(f1, abs(X1)/N1)
xlabel('Frequency (Hz)');
ylabel('Modulus')
axis([0 100 0 0.55])
subplot(3,1,3);
plot(f2, abs(X2)/N2)
xlabel('Frequency (Hz)');
ylabel('Modulus')
axis([0 50 0 0.55])