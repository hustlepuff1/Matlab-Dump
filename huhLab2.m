clc; close all;

test1_t = data(:,1);
test1_N = data(:,2);
test2_t = data(:,3);
test2_N = data(:,4);

figure('Name','B6-4 Test 1');
grid on
plot(test1_t,test1_N,'r');
title('B6-4 Test 1');
xlabel('Time(s)');
ylabel('Thrust(N)');

figure('Name','B6-4 Test 2');
grid on
plot(test2_t,test2_N,'b');
title('B6-4 Test 2');
xlabel('Time(s)');
ylabel('Thrust(N)');

figure('Name','B6-4 Test 1,2');
grid on
plot(test1_t,test1_N,'r');
hold on
plot(test2_t,test2_N,'b');
title('B6-4 Test 1, 2');
xlabel('Time(s)');
ylabel('Thrust(N)');
legend('Test 1','Test 2');