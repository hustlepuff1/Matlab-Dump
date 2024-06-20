clc; close all; clear all;
% Define parameters
A = 1;       % Decay constant
T = 0.01;  % Time increment
T1 = 10;     % Upper bound for time axis

% Generate time vector from 0 to T1 with increment T
t = 0:T:T1;

% Generate exponentially decaying signal
x = exp(-A * t);

% Define frequency axis for FFT
w = t;
w1 = -fliplr(w);
w2 = [w1 w];

% Compute Fourier Transform
ft = x * exp(-1i * t' * w) * T;
fft = fft(x);

% Compute magnitude and phase spectrum
a = abs(ft);
a1 = fliplr(a);
a2 = [a1 a];

p = angle(ft);
p1 = -fliplr(p);
p2 = [p1 p];

% Plotting the signal, magnitude spectrum, and phase spectrum
figure;

% Plot the exponential signal
subplot(3,1,1);
plot(t, x);
xlabel('Time');
ylabel('Amplitude');
grid on;
title('Exponential Signal');

% Plot the magnitude spectrum
subplot(3,1,2);
plot(w2, a2);
xlabel('Frequency');
ylabel('Magnitude');
grid on;
title('Magnitude Spectrum');

% Plot the phase spectrum
subplot(3,1,3);
plot(w2, p2);
xlabel('Frequency');
ylabel('Phase Angle (radians)');
grid on;
title('Phase Spectrum');
