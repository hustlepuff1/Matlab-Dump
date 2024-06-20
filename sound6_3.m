clc; clear all; close all;
f1=9; f2=14; fs=50; T=15.6;
t=0:1/fs:T-1/fs;
x=1*sin(2*pi*f1*t) + 0.001*sin(2*pi*f2*t); % Expression of the above equation
N=length(x);
whan=hann(N)'; % Generate the Hanning window
xh=x.*whan; % Apply the Hanning window
X=fft(x); % FFT of the original signal
Xh=fft(xh); % FFT of the windowed signal
f=fs*(0:N-1)/N; % Frequency vector

figure(1)
plot(f(1:N/2+1), 20*log10(abs(X(1:N/2+1)/fs/T)), 'r'); % Plot FFT of original signal
hold on
plot(f(1:N/2+1), 20*log10(sqrt(8/3)*abs(Xh(1:N/2+1)/fs/T)),'b'); % Plot FFT of windowed signal
axis([0 25 -180 0])
xlabel('Frequency (Hz)');
ylabel('Modulus (dB)')
legend({'Rectangular', 'Hann'});
hold off
