clc; close all;

of = cea(:,1);
ch4 = cea(:,2);
rp1 = cea(:,3);
h2 = cea(:,4);

figure('Name','ch4');
grid on
plot(of,ch4,'r');
title('Lox/LCH_{4}')
xlabel('Mixture(O/F) Ratio');
ylabel('Temperature [K]');

figure('Name','rp1');
grid on
plot(of,rp1,'b');
title('Lox/RP-1')
xlabel('Mixture(O/F) Ratio');
ylabel('Temperature [K]');

figure('Name','h2');
grid on
plot(of,h2,'k');
title('Lox/LH_{2}')
xlabel('Mixture(O/F) Ratio');
ylabel('Temperature [K]');

figure('Name','ch4');
grid on
hold on
plot(of,ch4,'r');
plot(of,rp1,'b');
plot(of,h2,'k');
title('Lox/LCH_{4}, RP-1, LH_{2}')
xlabel('Mixture(O/F) Ratio');
ylabel('Temperature [K]');
legend('LCH_{4}', 'RP-1', 'LH_{2}');