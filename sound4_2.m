clc; close all; clear all;

% Define parameters
T = 0.01;    % Time increment
T1 = 100;    % Upper bound for time axis

% Generate time vector from -T1 to T1 with increment T
t = -T1:T:T1-T;

% Generate rectangular pulse signal
x = rectpuls(t, 1);

% Compute Fourier Transform using explicit calculation
N = length(t);  % Length of time vector
Fs = 1/T;       % Sampling frequency
w = linspace(-Fs/2, Fs/2, N);  % Frequency axis centered at 0

ft = zeros(size(w));  % Initialize Fourier Transform
for k = 1:N
    ft(k) = sum(x .* exp(-1i * w(k) * t)) * T;
end


% Plotting the signal and magnitude spectrum
figure;

% Plot the rectangular pulse signal
subplot(2,1,1);
plot(t, x);
xlabel('Time');
ylabel('Amplitude');
xlim([-3, 3]);
ylim([-0.1, 1.1]);
grid on;
title('Rectangular Pulse Signal');

% Plot the magnitude spectrum
subplot(2,1,2);
plot(w, ft);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
ylim([-0.5, 1.2]);
grid on;
title('Magnitude Spectrum');
