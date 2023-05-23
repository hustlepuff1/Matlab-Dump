clc; close all;

moment = air5(:,1);
gauge1 = air5(:,2);
gauge2 = air5(:,3);
std1 = air5(:,4);
std2 = air5(:,5);

mdl1 = fitlm(gauge1, moment);
mdl2 = fitlm(gauge2, moment);

scatter(gauge1, moment, 'b');  % Scatter plot of gauge1 vs. moment in blue
hold on
scatter(gauge2, moment, 'r');  % Scatter plot of gauge2 vs. moment in red

% Plot fitted lines
x1 = min(gauge1):0.01:max(gauge1);
x2 = min(gauge2):0.01:max(gauge2);

y1 = predict(mdl1, x1');
plot(x1, y1, 'b', 'LineWidth', 2);  % Fitted line for mdl1 in blue

y2 = predict(mdl2, x2');
plot(x2, y2, 'r', 'LineWidth', 2);  % Fitted line for mdl2 in red

% Add error bars for standard deviation
errorbar(gauge1, moment, std1, 'horizontal', 'b.', 'LineWidth', 1);  % Error bars for gauge1
errorbar(gauge2, moment, std2, 'horizontal', 'r.', 'LineWidth', 1);  % Error bars for gauge2

legend('Gauge1', 'Gauge2', 'Model1', 'Model2', 'Std1', 'Std2');
xlabel('Voltage(V)');
ylabel('Bending moment(N-m)');
title('Moment - Voltage');

xlim([0, max([gauge1; gauge2])]);  % Set x-axis limits
ylim([0, max(moment)]);  % Set y-axis limits

hold off
