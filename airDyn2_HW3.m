clc; close all;

% 변수 지정
m1 = 1:0.1:10;
gam = 1.4;      % 공기 비열비
r = 287;        % 공기기체상수 R
Cp = (gam / (gam - 1)) * r;

% 식
m2 = sqrt((1 + ((gam - 1) / 2) * m1.^2) ./ (gam .* m1.^2 - (gam - 1) / 2));
rho2_rho1 = ((gam + 1) * m1.^2) ./ (2 + (gam - 1) * m1.^2);
p2_p1 = 1 + ((2 * gam) / (gam + 1)) * (m1.^2 - 1);
t2_t1 = p2_p1 ./ rho2_rho1;
p02_p01 = exp((-(Cp * log(t2_t1) - r * log(p2_p1))) / r);

% plot 왼쪽
yyaxis left;
grid on;
hold on;
plot(m1, m2, '-r');
plot(m1, p02_p01, '-g');
hold off;
ylabel('M_{2} and p_{0,2}/p_{0,1}')

% plot 오른쪽
yyaxis right;
hold on;
plot(m1, t2_t1, '-b');
plot(m1, p2_p1, '-k');
plot(m1, rho2_rho1, '-m');
hold off;
ylabel('T_{2}/T_{1}, p_{2}/p_{1}, and \rho_{2}/\rho_{1}');

% 그래프 꾸미기
ylim([0, 20]);
yticks(0:2:20);

legend('M_{2}', 'p_{0,2}/p_{0,1}', 'T_{2}/T_{1}', 'p_{2}/p_{1}', '\rho_{2}/\rho_{1}');
title('Normal shockwave properties. \gamma = 1.4');
xlabel('M_{1}');
% ylabel은 각자




