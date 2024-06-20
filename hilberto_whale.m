% Read the WAV file
[y, Fs] = audioread('whale.wav');

% Apply the Hilbert transform
y_hilbert = hilbert(y);

% Plot the original and transformed signals
t = (0:length(y)-1) / Fs; % Time vector
subplot(2,1,1);
plot(t, y);
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, abs(y_hilbert));
title('Absolute Value of Hilbert Transform');
xlabel('Time (s)');
ylabel('Amplitude');

% Optionally, you can play the original and transformed signals
% sound(y, Fs);
% sound(abs(y_hilbert), Fs);
