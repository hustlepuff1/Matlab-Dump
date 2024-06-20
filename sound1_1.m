%% 7음계의 monotone soundwave
clc; close all;

% 계이름 별 주파수
frequencies = [261.63, 293.66, 329.63, 349.23, 392.00, 440.00, 493.88, 523.25];

duration = 0.1; % 초
fs = 44100;   % Hz
t = linspace(0, duration, duration * fs);

% 사인파 영행렬 정의 후 만들기
note_sounds = zeros(length(frequencies), duration * fs);
for i = 1:length(frequencies)
    note_sounds(i, :) = sin(2 * pi * frequencies(i) * t);
end

 % 사인파 그리기
figure('Position', [100, 100, 600, 800]);
for j = 1:length(frequencies)
    subplot(length(frequencies), 1, j);
    plot(t, note_sounds(j, :)); % Plot the entire length of t
    title(sprintf('Note %d - Frequency %.2f Hz', j, frequencies(j)));
    xlabel('Time (s)');
    ylabel('Amp');
    ylim([-1.2, 1.2]);
    grid on;
end

set(gcf, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);

% 최대값 1, -1 로 제한
melody = [note_sounds(1,:) note_sounds(2,:) note_sounds(3,:) note_sounds(4,:) ...
          note_sounds(5,:) note_sounds(6,:) note_sounds(7,:) note_sounds(8,:)];
melody = melody / max(abs(melody));

% 소리 출력
sound(melody, fs);






