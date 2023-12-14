clc; close all;

mach = linspace(0, 10, 1001);
gam = 1.4;
rho = (1 + ((gam - 1) / 2) * mach.^2).^(-1 / (gam - 1));

plot(mach, rho, 'r.-');
grid on;

title('Density - Mach number relation');
xlabel('M');
ylabel('ρ/ρ_{0}');
