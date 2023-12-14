clc; close all;

alpha = -2:0.01:15;

figure('Name','C_{m} - \alpha', 'Position', [100, 100, 600, 800]);

no_It = -0.01 * alpha - 0.015;
plot(alpha, no_It, 'DisplayName', 'NO i_{t}');

hold on

for alpha_It = 2:1:13
    It = 0.75 + 0.6 * alpha_It;
    Cm = -0.01 * alpha + 0.02 * It - 0.015;
    plot(alpha, Cm, 'DisplayName', sprintf('i_{t}[%d]', alpha_It));
end

legend('Location', 'northeastoutside');
xlim([-2, 14]);
ylim([-0.25, 0.25]);
grid on;
title('C_{m} - \alpha');
xlabel('\alpha');
ylabel('C_{m}');
grid on;

It_table = array2table(It_values', 'VariableNames', {'It Values'});
disp(It_table);

