clc; close all;

theta = data;
speed5 = data1;
speed30 = dataS2;

theta2 = linspace(0,360,36);
radian = deg2rad(theta);

speed51 = abs(speed5/9.81);
speed301 = abs(speed30/9.81);

%y = 1/2*(1.23*5^2)*(1-4*sin(radian).^2);

%polarplot(y, 'b.-');

% hold on;
% r = 10;
% theta_circle = linspace(0, 2*pi, 100);
% rho_circle = ones(size(theta_circle)) * r;
% polarplot(theta_circle, rho_circle, 'b', 'LineWidth', 2);

subplot(2,2,1)
plot(theta,speed5,"r.-");
subplot(2,2,2);
polarplot(radian,speed51,"r.-");
subplot(2,2,3);
plot(theta,speed30,"b.-");
subplot(2,2,4);
polarplot(radian,speed301,"b.-");


