clear; clc; close;

[X,Y] = meshgrid(-3:0.5:3);
a = cos(X) .* cosh(Y);
vx = -sin(X) .* cosh(Y);
vy = cos(X) .* sinh(Y);

hold on

h = quiver(X, Y, vx, vy);
set(h, 'color', 'g', 'linewidth', 2);
contour(X, Y, a);

grid on

hold off