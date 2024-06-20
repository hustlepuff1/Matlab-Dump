clc; clear all; close all;

fs = 100; 
t = 0:1/fs:5-1/fs;
A = 200; 
zeta = 0.01; 
wn = 2*pi*20;
wd = sqrt(1-zeta^2)*wn;
x = (A/wd)*exp(-zeta*wn*t).*sin(wd*t); % Expression of the time signal
N = length(x);
whan = hann(N)'; % Generate the Hanning window
xh = x.*whan; % Apply the Hanning window
X = fft(x); % FFT of the original signal
Xh = fft(xh); % FFT of the windowed signal
f = fs*(0:N-1)/N; % Frequency vector

H = A./(wn^2 - (2*pi*f).^2 + 1i*2*zeta*wn*(2*pi*f)); % Transfer function

subplot(2,1,1)
plot(f(1:N/2+1), 20*log10(abs(X(1:N/2+1)/fs))); % Plot FFT of original signal
hold on
plot(f(1:N/2+1), 20*log10(sqrt(8/3)*abs(Xh(1:N/2+1)/fs)), 'r'); % Plot FFT of windowed signal
plot(f(1:N/2+1), 20*log10(abs(H(1:N/2+1))), 'k'); % Plot transfer function
axis([0 50 -150 0])
xlabel('Frequency (Hz)'); 
ylabel('Modulus (dB)')
hold off

subplot(2,1,2)
plot(f(1:N/2+1), abs(X(1:N/2+1)/fs)); 
hold on
plot(f(1:N/2+1), sqrt(8/3)*abs(Xh(1:N/2+1)/fs), 'r');
plot(f(1:N/2+1), abs(H(1:N/2+1)), 'k');
axis([0 50 0 0.7])
xlabel('Frequency (Hz)');
ylabel('Modulus (linear scale)')
legend({'Rectangular window', 'Hann window', 'True magnitude spectrum'});
hold off
