clc; close all; clear all;

% Define constants
f1 = 10; 
f2 = 20; 
f3 = 21; 
fs = 60;
T_values_1 = [0.6, 0.8, 1.0, 1.5]; % First set of T values
T_values_2 = [2, 2.5, 3, 4]; % Second set of T values

% Define number of subplots for each figure
num_plots_1 = length(T_values_1);
num_plots_2 = length(T_values_2);

% Determine the grid size for subplots
num_rows = 4; % Set the number of rows to 4
num_cols = 1; % Set the number of columns to 1

% Create the first figure for the first set of T values
figure;
% Loop through each T value in the first set
for idx = 1:num_plots_1
    T = T_values_1(idx);
    
    % Compute time vector
    t = 0:1/fs:T-1/fs;
    
    % Generate the signal
    x = 2*sin(2*pi*f1*t) + 2*sin(2*pi*f2*t) + 2*sin(2*pi*f3*t);
    
    % Length of the signal
    N = length(x);
    
    % Compute FFT for original signal (without window)
    X = fft(x);
    
    % Frequency vector for original signal
    f = fs * (0:N-1) / N;
    
    % Zero-padding and FFT for original signal
    Xz = fft([x zeros(1,2000-N)]);
    Nz = length(Xz);
    fz = fs * (0:Nz-1) / Nz;
    
    % Apply Hanning window
    whan = hanning(N);
    x_win = x .* whan';
    
    % Compute FFT for windowed signal
    X_win = fft(x_win);
    
    % Zero-padding and FFT for windowed signal
    Xz_win = fft([x_win zeros(1,2000-N)]);
    Nz_win = length(Xz_win);
    fz_win = fs * (0:Nz_win-1) / Nz_win;
    
    % Create a subplot
    subplot(num_rows, num_cols, idx);
    
    % Plot results for original signal in red
    h1 = stem(f(1:N/2+1), abs(X(1:N/2+1)/fs/T), 'r:', 'DisplayName', 'Rectangular');
    hold on;
    h2 = plot(fz(1:Nz/2+1), abs(Xz(1:Nz/2+1)/fs/T), 'r', 'DisplayName', 'Rectangular (Zero-padded)');
    
    % Plot results for windowed signal in blue
    h3 = stem(f(1:N/2+1), sqrt(8/3)*abs(X_win(1:N/2+1)/fs/T), 'b:', 'DisplayName', 'Hann');
    h4 = plot(fz_win(1:Nz_win/2+1), sqrt(8/3)*abs(Xz_win(1:Nz_win/2+1)/fs/T), 'b', 'DisplayName', 'Hann (Zero-padded)');
    hold off;
    
    % Add labels, axis limits, and grid
    title(sprintf('T = %.1f', T));
    xlabel('Frequency (Hz)');
    ylabel('Modulus');
    axis([0 30 0 1.2]);
    grid on;
    
    % Store legend entries for the first subplot only
    if idx == 1
        legend_entries = [h1, h2, h3, h4];
    end
end

% Add a single legend for the entire first figure
legend(legend_entries, {'Rectangular', 'Rectangular (Zero-padded)', 'Hann', 'Hann (Zero-padded)'}, 'Location', 'bestoutside');

% Create the second figure for the second set of T values
figure;
% Loop through each T value in the second set
for idx = 1:num_plots_2
    T = T_values_2(idx);
    
    % Compute time vector
    t = 0:1/fs:T-1/fs;
    
    % Generate the signal
    x = 2*sin(2*pi*f1*t) + 2*sin(2*pi*f2*t) + 2*sin(2*pi*f3*t);
    
    % Length of the signal
    N = length(x);
    
    % Compute FFT for original signal (without window)
    X = fft(x);
    
    % Frequency vector for original signal
    f = fs * (0:N-1) / N;
    
    % Zero-padding and FFT for original signal
    Xz = fft([x zeros(1,2000-N)]);
    Nz = length(Xz);
    fz = fs * (0:Nz-1) / Nz;
    
    % Apply Hanning window
    whan = hanning(N);
    x_win = x .* whan';
    
    % Compute FFT for windowed signal
    X_win = fft(x_win);
    
    % Zero-padding and FFT for windowed signal
    Xz_win = fft([x_win zeros(1,2000-N)]);
    Nz_win = length(Xz_win);
    fz_win = fs * (0:Nz_win-1) / Nz_win;
    
    % Create a subplot
    subplot(num_rows, num_cols, idx);
    
    % Plot results for original signal in red
    h1 = stem(f(1:N/2+1), abs(X(1:N/2+1)/fs/T), 'r:', 'DisplayName', 'Rectangular');
    hold on;
    h2 = plot(fz(1:Nz/2+1), abs(Xz(1:Nz/2+1)/fs/T), 'r', 'DisplayName', 'Rectangular (Zero-padded)');
    
    % Plot results for windowed signal in blue
    h3 = stem(f(1:N/2+1), sqrt(8/3)*abs(X_win(1:N/2+1)/fs/T), 'b:', 'DisplayName', 'Hann');
    h4 = plot(fz_win(1:Nz_win/2+1), sqrt(8/3)*abs(Xz_win(1:Nz_win/2+1)/fs/T), 'b', 'DisplayName', 'Hann (Zero-padded)');
    hold off;
    
    % Add labels, axis limits, and grid
    title(sprintf('T = %.1f', T));
    xlabel('Frequency (Hz)');
    ylabel('Modulus');
    axis([0 30 0 1.2]);
    grid on;
    
    % Store legend entries for the first subplot only
    if idx == 1
        legend_entries = [h1, h2, h3, h4];
    end
end

% Add a single legend for the entire second figure
legend(legend_entries, {'Rectangular', 'Rectangular (Zero-padded)', 'Hann', 'Hann (Zero-padded)'}, 'Location', 'bestoutside');
