clc; close all;
data = readmatrix('cl.csv');
x1 = data(:,1);
y1 = data(:,2);
x2 = data(:,4);
y2 = data(:,5);

plot(x1,y1, "Color", 'g');
hold on
plot(x2,y2, "Color","#77AC30");
xlabel('Alpha');
ylabel('CL');
title('CL - Alpha Graph (Inviscid)');
grid on

% find and mark max point of the green line
[max_y1, max_idx1] = max(y1);
max_x1 = x1(max_idx1);
plot(max_x1, max_y1, 'ro', 'MarkerSize', 10, 'LineWidth', 2)

% find and mark max point of the red line
[max_y2, max_idx2] = max(y2);
max_x2 = x2(max_idx2);
plot(max_x2, max_y2, 'bo', 'MarkerSize', 10, 'LineWidth', 2)
legend('Re: 1.7e3','Re: 3.55e5', 'Stall at Re: 1.7e3','Stall at Re: 3.55e5','Location','northwest');
