clear all; close all; clc
% TF of elements
sys_coupler = 10*tf(conv([1 0.1],[1 0.5]),conv([1 0],[1 5])) % coupler
sys_pitch = 46.5*tf([1 0.585],conv([1 6 0.5*5.5],[1 5.4 11.4])) % pitch
sys_gamma = -0.0276*tf(conv([1 4.85],[1 -4.35]),[1 0.585]) % gamma/theta
u0 = 280;
sys_d = tf(u0, 57.3*[1 0]) % d/gamma
sys_Gamma = 57.3 % Gamma/d
% root locus
a = 10*46.5*0.0276*u0
G = sys_coupler*sys_pitch*sys_gamma*sys_d*sys_Gamma/a % forward
H = 1 % feedback
GH = G*H % for root locus (1+kGH=0)
rlocus(GH); grid;
% impulse response
k = 10;
cltf = feedback(k*G,H);
figure; impulse(cltf)