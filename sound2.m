clc; close all;

A1 = 1;
A2 = 1;
Theta1 = 0;
Theta2 = 0;
f1_1 = 1.4;
f1_2 = 1.5;
f2_1 = sqrt(2);
f2_2 = 1.5;

t = 0:0.01:30;

x1 = A1*sin(2*pi*f1_1*t + Theta1) + A2*sin(2*pi*f1_2*t + Theta2);
x2 = A1*sin(2*pi*f2_1*t + Theta1) + A2*sin(2*pi*f2_2*t + Theta2);

X1 = fft(x1);
f1 = linspace(0, 1/(2*(t(2)-t(1))), floor(length(X1)/2));
magnitude_X1 = abs(X1(1:length(X1)/2));

X2 = fft(x2);
f2 = linspace(0, 1/(2*(t(2)-t(1))), floor(length(X2)/2));
magnitude_X2 = abs(X2(1:length(X2)/2));

figure(1)
subplot(2, 1, 1);
plot(t, x1, 'b');
xlabel('Seconds')
ylabel('x(t)')
title('Problem 1 (1.4 Hz + 1.5 Hz)')
grid on

subplot(2, 1, 2);
plot(f1, magnitude_X1, 'b');
title('Problem 1 - FFT (1.4 Hz + 1.5 Hz)')
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 3]);
grid on

figure(2);
subplot(2, 1, 1);
plot(t, x2, 'r');
xlabel('Seconds')
ylabel('x(t)')
title('Problem 2 (\surd{2} Hz + 1.5 Hz)')
grid on

subplot(2, 1, 2);
plot(f2, magnitude_X2, 'r');
title('Problem 2 - FFT (\surd{2} Hz + 1.5 Hz)')
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([0 3]);
grid on

figure(3);
hold on
plot(f1, magnitude_X1, 'b');
plot(f2, magnitude_X2, 'r');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Problem 1 (1.4 Hz + 1.5 Hz) + Problem 2 (\surd{2} Hz + 1.5 Hz)')
xlim([0 3]);
legend('1.4 Hz + 1.5 Hz', '\surd{2} Hz + 1.5 Hz', 'Location', 'NorthEast');
grid on

sound(x1);
sound(x2);
