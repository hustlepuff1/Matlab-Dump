clear all; close all; clc;

m = 0.111;
R = 0.015;
g = -9.8;
L = 1.0;
d = 0.03;
J = 9.99e-6;

s = tf('s');

P_ball = -(m*g*d)/(L*(J/(R^2)+m))*1/s^2;

Kp = 1;
C = pid(Kp);
sys_cl = feedback(C*P_ball, 1);

Kp1 = 1; Kd1 = 1;
Kp2 = 10; Kd2 = 10;
Kp3 = 10; Kd3 = 20;
Kp4 = 15; Kd4 = 40;

C1 = pid(Kp1, 0, Kd1);
C2 = pid(Kp2, 0, Kd2);
C3 = pid(Kp3, 0, Kd3);
C4 = pid(Kp4, 0, Kd4);

sys_cl1 = feedback(C1*P_ball, 1);
sys_cl2 = feedback(C2*P_ball, 1);
sys_cl3 = feedback(C3*P_ball, 1);
sys_cl4 = feedback(C4*P_ball, 1);

figure(1)
pzmap(P_ball)
figure(2)
step(0.25*P_ball)

figure(3)
step(0.25*sys_cl)
axis([0 70 0 0.5])

figure(4)
for i = 1:4
    subplot(2, 2, i);
    switch i
        case 1
            sys_cl_i = sys_cl1;
            title('Controller 1');
        case 2
            sys_cl_i = sys_cl2;
            title('Controller 2');
        case 3
            sys_cl_i = sys_cl3;
            title('Controller 3');
        case 4
            sys_cl_i = sys_cl4;
            title('Controller 4');
    end
    step(0.25 * sys_cl_i);
    axis([0 70 0 0.5]);
    legend('Location', 'northeast');

    step_info = stepinfo(sys_cl_i);
    fprintf('Controller %d Step Response Info:\n', i);
    disp(step_info);
end




