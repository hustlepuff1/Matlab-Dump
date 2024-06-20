% Clear all previous data
clear all

% Define parameters
A = 2; 
p = 1; 
Tp = 1 / p; 
fs = 10 / Tp;

% Define truncation periods
T1 = 1.5 * Tp; 
T2 = 3.5 * Tp;

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

% Perform 5000-point DFT by adding zeros at the end of each sequence
X1z = fft([x1 zeros(1,5000-N1)]); % zero padding
X2z = fft([x2 zeros(1,5000-N2)]); % zero padding

% Calculate the length of the zero-padded DFT
Nz = length(X1z);

% Calculate the new frequency variable
fz = fs * (0:Nz-1) / Nz;

% Plot the results (modulus) of 15-point DFT and DFT with zero padding
subplot(2,1,1);
stem(f1, abs(X1) / fs / T1, 'fill');
hold on
plot(fz, abs(X1z) / fs / T1, 'r'); 
hold off
xlabel('Frequency (Hz)');
ylabel('Modulus (scaled)')
axis([0 10 0 1.02])

% Plot the results (modulus) of 35-point DFT and DFT with zero padding
subplot(2,1,2);
stem(f2, abs(X2) / fs / T2, 'fill'); 
hold on
plot(fz, abs(X2z) / fs / T2, 'r'); 
hold off
xlabel('Frequency (Hz)');
ylabel('Modulus (scaled)')
axis([0 10 0 1.02])
