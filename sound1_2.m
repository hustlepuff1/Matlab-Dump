%% 음원파일 FFT 하기
% 음원파일 가져오기
[y, Fs] = audioread('wow.wav');
sound(y, Fs)

% 길이 알기
nfft = length(y);

% Compute FFT
Y = fft(y, nfft);
Y = Y(1:round(nfft/2)); % 양수값만 추출

% 시간 백터
t = (0:length(y)-1) / Fs;
% 주파수 백터 계산
f = Fs*(0:round(nfft/2)-1)/nfft;

% Plot
subplot(2, 1, 1);
plot(t, y);
xlabel('Time (s)');
ylabel('Amplitude');
title('Time domain');

subplot(2, 1, 2);
plot(f, abs(Y));
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FFT Result');
