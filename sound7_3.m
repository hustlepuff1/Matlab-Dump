clc; clear all; close all;

fs = 500;
T = 10;
t = 0:1/fs:T-1/fs;
p1 = 10;
p2 = 40;

x = sin(2*pi*p1*t) + sin(2*pi*p2*t);

x1 = resample(x, 100, fs);
x2 = resample(x, 50, fs);

N = length(x);
N1 = length(x1);
N2 = length(x2);

X = fft(x);
X1 = fft(x1);
X2 = fft(x2);

f = fs*(0:N-1)/N;
f1 = 100*(0:N1-1)/N1;
f2 = 50*(0:N2-1)/N2;

subplot(3,1,1);
plot(f, abs(X)/N)
xlabel('Frequency (Hz)');
ylabel('Modulus')
axis([0 500 0 0.55])
title('Frequency Spectrum of Original Signal')
subplot(3,1,2);
plot(f1, abs(X1)/N1)
xlabel('Frequency (Hz)');
ylabel('Modulus')
axis([0 100 0 0.55])
title('Frequency Spectrum of Resampled Signal at 100 Hz')
subplot(3,1,3);
plot(f2, abs(X2)/N2)
xlabel('Frequency (Hz)');
ylabel('Modulus')
axis([0 50 0 0.55])
title('Frequency Spectrum of Resampled Signal at 50 Hz')