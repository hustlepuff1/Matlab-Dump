clc; clear; close all;

t = 0:0.1:10; % each spaced 0.1s apart
alpha = 1;

x = exp(-alpha*t); % x is now a list of raw data


L = length(x); % number of samples (samples)
Fs = mean(diff(t)); % average sampling frequency between...
                    % ...my raw data (samples/sec)
f = (-L/2:L/2 -1) * Fs/L; % my x-axis =>
                    % frequency (Hz, 1/s) samples/sec / samples

% plotting my raw data
figure(1);
plot(t, x);
xlabel('Time (s)');
title('Exponential');
grid on

% trying a FFT
X = fft(x); % large X is just a complex number!
            % if we plot this, it gives us garbage! :(
            % angle and amplitude
X_norm = 1/L * X; % lets normalize by divide the number of terms

% Amplitude and phase spectrum
figure(2)
subplot(1,2,1)
plot(f, fftshift(abs(X_norm)))
xlabel('frequency (Hz)')
ylabel('Amplitude')
grid on

subplot(1,2,2)
fft_rad = angle(X_norm);
fft_ang = rad2deg(fft_rad);
plot(f, fftshift(fft_rad))
xlabel('frequency (Hz)')
ylabel('Phase Angle')
grid on

% Make a table of values (coedds, freq, amp, ang)
table1 = table(X_norm',f',abs(X_norm'),fft_ang');
table1.Properties.VariableNames = {'FFT Coeffs','Frequency','Amplitude','Phase'};
disp(table1)


