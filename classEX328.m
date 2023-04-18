num = 30*[1 -6];
den = [1 4 13 0];
sys = tf(num, den)*(-1);

t = 0:0.01:5;
impulse(sys, t);
t = 0:0.01:0.7;
figure; step(sys, t);