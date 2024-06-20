% Clear all previous data
clear all

% Define parameters
A = 2; 
p = 1; 
Tp = 1 / p; 
fs = 20 / Tp;

% Define truncation periods
T1 = 1 * Tp; 
T2 = 1.5 * Tp;

% Define time variables
t1 = [0:1/fs:T1-1/fs]; 
t2 = [0:1/fs:T2-1/fs];

% Generate sampled and truncated signals
x1 = A * cos(2 * pi * p * t1);
x2 = A * cos(2 * pi * p * t2);

% Perform the DFT of each signal
X1 = fft(x1);
X2 = fft(x2);

% Calculate the lengths of each signal
N1 = length(x1); 
N2 = length(x2);

% Calculate the frequency variables
f1 = fs * (0:N1-1) / N1; 
f2 = fs * (0:N2-1) / N2;

% Plot the results (modulus) of 20-point DFT
subplot(2,1,1);
stem(f1, abs(X1) / fs / T1, 'fill')
xlabel('Frequency (Hz)')
ylabel('Modulus (scaled)')
axis([0 20 0 1])

% Plot the results (modulus) of 30-point DFT
subplot(2,1,2);
stem(f2, abs(X2) / fs / T2, 'fill')
xlabel('Frequency (Hz)')
ylabel('Modulus (scaled)')
axis([0 20 0 1])
