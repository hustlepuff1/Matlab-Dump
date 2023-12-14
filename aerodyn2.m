clc; close all;

alpha = -2:0.01:15;

figure('Name','C_{m} - \alpha with \delta_{ele}', 'Position', [100, 100, 600, 800]);

del_values = zeros(1, numel(2:1:13));

for alpha_del = 2:1:13
    It = 3.75;
    del = -0.95 * alpha_del + 4.75;
    del_values(alpha_del - 1) = del;
    Cm = -0.01 * alpha + 0.02 * It - 0.015 - 0.01579 * del;
    plot(alpha, Cm, 'DisplayName', sprintf('\\delta_{ele}[%.2f]', del));
    hold on
end

legend('Location', 'northeastoutside');
xlim([-2, 14]);
ylim([-0.25, 0.25]);
grid on;
title('C_{m} - \alpha with \delta_{ele}');
xlabel('\alpha');
ylabel('C_{m}');
grid on;

del_table = array2table(del_values', 'VariableNames', {'delta_E_Values'});
disp(del_table);
