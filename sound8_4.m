% Clear all previous data
clear all

% Define the parameter
a = 0.3;

% Define the variables for the sequences
n1 = 0:8; % 9-point sequence
n2 = 0:9; % 10-point sequence

% Create the sequences
x1 = a.^n1; 
x2 = a.^n2;

% Perform the DFT of each sequence
X1 = fft(x1);
X2 = fft(x2);

% Plot the real part of the DFT of the first sequence x1
figure(1)
subplot(2,2,1);
stem(n1, real(X1), 'fill')
axis([-0.5 8.5 0 1.6])
xlabel('\itk');
ylabel('Re[\itX\rm(\itk\rm)]')

% Plot the imaginary part of the DFT of the first sequence x1
subplot(2,2,2);
stem(n1, imag(X1), 'fill')
axis([-0.5 8.5 -0.4 0.4])
xlabel('\itk');
ylabel('Im[\itX\rm(\itk\rm)]')

% Plot the modulus of the DFT of the first sequence x1
subplot(2,2,3);
stem(n1, abs(X1), 'fill')
axis([-0.5 8.5 0 1.6])
xlabel('\itk');
ylabel('|\itX\rm(\itk\rm)|')

% Plot the phase of the DFT of the first sequence x1
subplot(2,2,4);
stem(n1, angle(X1), 'fill')
axis([-0.5 8.5 -0.4 0.4])
xlabel('\itk');
ylabel('arg\itX\rm(\itk\rm)')

% Plot the real part of the DFT of the second sequence x2
figure(2)
subplot(2,2,1);
stem(n2, real(X2), 'fill')
axis([-0.5 9.5 0 1.6])
xlabel('\itk');
ylabel('Re[\itX\rm(\itk\rm)]')

% Plot the imaginary part of the DFT of the second sequence x2
subplot(2,2,2);
stem(n2, imag(X2), 'fill')
axis([-0.5 9.5 -0.4 0.4])
xlabel('\itk');
ylabel('Im[\itX\rm(\itk\rm)]')

% Plot the modulus of the DFT of the second sequence x2
subplot(2,2,3);
stem(n2, abs(X2), 'fill')
axis([-0.5 9.5 0 1.6])
xlabel('\itk');
ylabel('|\itX\rm(\itk\rm)|')

% Plot the phase of the DFT of the second sequence x2
subplot(2,2,4);
stem(n2, angle(X2), 'fill')
axis([-0.5 9.5 -0.4 0.4])
xlabel('\itk');
ylabel('arg\itX\rm(\itk\rm)')
