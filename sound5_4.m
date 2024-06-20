clc; clear all; close all;

% 음원파일 가져오기
[y, Fs] = audioread('sound.mp3');
sound(y, Fs)

cat_hilbert = hilbert(y);

% 원음 plot
t = (0:length(y)-1) / Fs;
subplot(3,1,1);
plot(t, y);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Audio Signal');
grid on;

% Hilbert transform
subplot(3,1,2);
hold on
plot(t, y, 'r');
plot(t, imag(cat_hilbert), 'b');
xlabel('Time (s)');
ylabel('Amplitude');
title('Analytic Signal (Hilbert Transform)');
grid on;

subplot(3,1,3);
hold on
plot(t, y, 'r');
plot(t, imag(cat_hilbert), 'b');
xlabel('Time (s)');
ylabel('Amplitude');
title('Analytic Signal (Hilbert Transform)');
axis([0.23 0.25 -1 1])
grid on;

sgtitle('Original Audio Signal and Analytic Signal (Hilbert Transform)');

