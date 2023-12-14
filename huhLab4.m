clc; close all;

t = data(:,1);
exp1 = data(:,2);
exp2 = data(:,3);
exp3 = data(:,4);

figure('Name','Thrust');
grid on
hold on
plot(t,exp1,'r');
plot(t,exp2,'b');
plot(t,exp3,'K');
xlabel('Time, s');
ylabel('Thrust, N');
