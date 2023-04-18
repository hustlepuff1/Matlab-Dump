clear; clc; close;

[X,Y,Z] = meshgrid(-10:2:10);
vx = exp(-Z .^ 2);
vy = exp(-X .^ 2);
vz = exp(-Y .^ 2);

h = quiver3(X, Y, Z, vx, vy, vz);
axis image
grid on
set(h, 'color', 'r', 'linewidth', 1);
hold on
[u,v,w] = curl(vx, vy, vz);
k = quiver3(X, Y, Z, u, v, w);
set(h, 'color', 'b', 'linewidth', 1);
hold off