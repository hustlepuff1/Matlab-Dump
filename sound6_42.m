clc; clear all; close all;
fs=100; t=0:1/fs:5-1/fs;
A=200; B=0.001*A; zeta1=0.01; zeta2=0.01;
wn1=2*pi*20; wd1=sqrt(1-zeta1^2)*wn1;
wn2=2*pi*30; wd2=sqrt(1-zeta2^2)*wn2;
x=(A/wd1)*exp(-zeta1*wn1*t).*sin(wd1*t) + (B/wd2)*exp(-zeta2*wn2*t).*sin(wd2*t);
N=length(x);
whan=hanning(N); xh=x.*whan';
X=fft(x); Xh=fft(xh);
f=fs*(0:N-1)/N;
H=A./(wn1^2-(2*pi*f).^2+i*2*zeta1*wn1*(2*pi*f)) + B./(wn2^2-(2*pi*f).^2+i*2*zeta2*wn2*(2*pi*f));

subplot(2,1,1)
plot(f(1:N/2+1), 20*log10(abs(X(1:N/2+1)/fs)), 'r');
hold on
plot(f(1:N/2+1), 20*log10(abs(H(1:N/2+1))), 'b')
axis([0 50 -60 0])
xlabel('Frequency (Hz)'); ylabel('Modulus (dB)')
legend({'Rectangular window', 'True magnitude spectrum'});
hold off

subplot(2,1,2)
plot(f(1:N/2+1), 20*log10(sqrt(8/3)*abs(Xh(1:N/2+1)/fs)), 'r')
hold on
plot(f(1:N/2+1), 20*log10(abs(H(1:N/2+1))), 'b')
axis([0 50 -160 0])
xlabel('Frequency (Hz)'); ylabel('Modulus (dB)')
legend({'Rectangular window', 'Hann window'});
hold off