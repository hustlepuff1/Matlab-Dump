clc; close all;

%1.2 문제
num = [1, 20];
den = [1, 8, 25];
Gtf = tf(num, den);
[z, p, k] = tf2zp(num, den);
fprintf('k = %g\n', k); 
fprintf('z = ');
fprintf('%g ', z);
fprintf('\n');
fprintf('p = ');
disp(p);
fprintf('\n');

%1.4 문제
t = 0:0.1:10;
% 단위계단
subplot(2,1,1);
step(Gtf, t);
% 임펄스
subplot(2,1,2);
impulse(Gtf, t);

%1.5 문제
syms s t;
Ys = (s+20)/(s^2+8*s+25);
% 단위계단
YsStep = (s+20)/(s*(s^2+8*s+25));
ytStep = ilaplace(YsUs);
disp(ytStep);
limStep = limit(ytStep, t, inf);
disp(ytStep);
disp(limStep);
% 임펄스
YsImp = (s+20)/(s^2+8*s+25);
ytImp = ilaplace(YsImp);
limImp = limit(ytImp, t, inf);
disp(ytImp);
disp(limImp);

%2.2 문제
% G1 정의
num1 = 25;
den1 = [1, 5, 25];
G1 = tf(num1 , den1);
% G2 정의
num2 = 1;
den2 = [1, 1];
G2 = tf(num2, den2);
% G3 Backward Inner Loop
num3 = 1;
den3 = 1;
G3 = tf(num3, den3);
% 최종 연산
sysInner = feedback(G1, G3);
sysA = parallel(sysInner, G2);