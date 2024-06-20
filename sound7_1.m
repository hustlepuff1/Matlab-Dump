clc; clear all; close all;

fs=100; T=10;
t=0:1/fs:T-1/fs;
p1=20; p2=80; p3=120;
x1=sin(2*pi*p1*t);
x2=sin(2*pi*p2*t);
x3=sin(2*pi*p3*t);
N=length(t);
X1=fft(x1); X2=fft(x2);
X3=fft(x3);
f=fs*(0:N-1)/N;
subplot(3,1,1); plot(f, abs(X1)/fs/T)
xlabel('Frequency (Hz)');
ylabel('Modulus')
axis([0 100 0 0.55])
subplot(3,1,2); plot(f, abs(X2)/fs/T)
xlabel('Frequency (Hz)');
ylabel('Modulus')
axis([0 100 0 0.55])
subplot(3,1,3); plot(f, abs(X3)/fs/T)
xlabel('Frequency (Hz)');
ylabel('Modulus')
axis([0 100 0 0.55])