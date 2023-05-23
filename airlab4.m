clc; close all;

height = data(:,1);
vel5 = data(:,10);
vel20 = data(:,11);
err5 = data(:,8);
err20 = data(:,9);

figure('name', 'Wake for 5deg');
plot(vel5, height, 'b');
errorbar(vel5, height, err5, 'both', 'b');
title('Wake at 5deg')
xlabel('U/U∞')
ylabel('Height(mm)')

figure('name', 'Wake for 20deg');
plot(vel20, height, 'r');
errorbar(vel20, height, err20, 'both', 'r');
title('Wake at 20deg')
xlabel('U/U∞')
ylabel('Height(mm)')