clc; clear all; close all;

% Define constants
f1 = 10; 
f2 = 20; 
f3 = 21;
fs = 60;
T_values = [0.6, 0.8, 1.0, 1.5, 2, 2.5, 3, 4]; % Different values of T

% Number of subplots
num_plots = length(T_values);

% Determine the grid size for subplots
num_rows = 4; % Set the number of rows to 4
num_cols = 2; % Set the number of columns to 2

% Create a new figure
figure;

% Loop through each T value
for idx = 1:num_plots
    T = T_values(idx);
    
    % Compute time vector
    t = 0:1/fs:T-1/fs;
    
    % Generate the signal
    x = 2*sin(2*pi*f1*t) + 2*sin(2*pi*f2*t) + 2*sin(2*pi*f3*t);
    
    % Length of the signal
    N = length(x);
    
    % Compute FFT
    X = fft(x);
    
    % Frequency vector
    f = fs * (0:N-1) / N;
    
    % Zero-padding and FFT
    Xz = fft([x zeros(1,2000-N)]);
    Nz = length(Xz);
    fz = fs * (0:Nz-1) / Nz;
    
    % Create a subplot
    subplot(num_rows, num_cols, idx);
    
    % Plot results
    stem(f(1:N/2+1), abs(X(1:N/2+1)/fs/T), 'r:');
    hold on;
    plot(fz(1:Nz/2+1), abs(Xz(1:Nz/2+1)/fs/T));
    hold off;
    
    % Add labels, axis limits, and grid
    title(sprintf('T = %.1f', T));
    xlabel('Frequency (Hz)');
    ylabel('Modulus');
    axis([0 30 0 1.2]);
    grid on;
end
