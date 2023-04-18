clear; clc; close;

[X,Y] = meshgrid(-3:0.5:3);
vx = X ./ (X .^ 2 + Y .^ 2);
vy = Y ./ (X .^ 2 + Y .^ 2);

hold on

h = quiver(X, Y, vx, vy);
set(h, 'color', 'r', 'linewidth', 2);

grid on

hold off