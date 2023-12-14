clc; close all;

alpha = data(:,1);
cl = data(:,2);
cd = data(:,3);
errcL = data(:,4);
errcD = data(:,5);
errcL2 = data(:,6);
errcL_lvm = data(:,7);

figure('name', 'cl-alpha');
plot(alpha, cl, 'b');
errorbar(alpha, cl, errcL_lvm, 'r');
title('Cl-Alpha')
xlabel('Alpha')
ylabel('Cl')

% figure('name', 'Wake for 20deg');
% plot(vel20, height, 'r');
% errorbar(vel20, height, err20, 'both', 'r');
% title('Wake at 20deg')
% xlabel('U/U∞')
% ylabel('Height(mm)')