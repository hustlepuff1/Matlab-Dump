clc; close all;

% 변수 설정, Import from plotting.xlsx
ang = 0:10:360;
angp = ang';
% c51 = Cp5;
% c301 = Cp30;
% cI = CpIdeal;
% c52 = abs(Cp5);
% c302 = abs(Cp30);
pa5 = abs(Pa5);
pa30 = abs(Pa30);

ang = 0:10:360;
rad = deg2rad(ang);


% --------------------------이상적인 압력------------------------------- %
figure(1);
y = 1/2*(1.23*5^2)*(1-4*sin(rad).^2);
polarplot(rad, y);
hold on
for i = 1:length(rad)
    polarplot([rad(i) rad(i)], [0 y(i)], 'k-');
end
% theta = linspace(0,2*pi,100);
% r = ones(1,length(theta))*0.05;
% polarplot(theta, r, 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'MarkerSize', 15);
hold off
title('Pressure Distribution at Ideal');

% ---------------------------차압 at 5m/s------------------------------ %
figure(2);
polarplot(rad, pa5, 'r');
hold on
for i = 1:length(rad)
    polarplot([rad(i) rad(i)], [0 pa5(i)], 'k-');
end
%theta = linspace(0,2*pi,100);
%r = ones(1,length(theta))*0.05;
%polarplot(theta, r, 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'MarkerSize', 15);
hold off
title('Pressure Distribution at 5m/s');

% ---------------------------차압 at 30m/s------------------------------ %
figure(3);
polarplot(rad, pa30, 'b');
hold on
for i = 1:length(rad)
    polarplot([rad(i) rad(i)], [0 pa30(i)], 'k-');
end
%theta = linspace(0,2*pi,100);
%r = ones(1,length(theta))*0.05;
%polarplot(theta, r, 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'MarkerSize', 15);
hold off
title('Pressure Distribution at 30m/s');



% -----------------------Plotting x,y Graph---------------------------- %
% figure(4);
% plot(angp, c51, angp, c301, angp, cI);
% xlim([0, 360]);
% xticks([0 90 180 270 360]);
% title('Cp, Angle Graph');
% xlabel('Angle');
% ylabel('Cp');
% legend('5 m/s', '30 m/s', 'Ideal');

% grid on

%%---------------------------------------------------------------------%%
