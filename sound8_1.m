clear all;
x=[1 1 1 1 1 0 0 0 0];
h=[8 7 6 5 4 3 2 1 0 0 0];
nx=[0:length(x)-1];
nh=[0:length(h)-1];
y1=conv(h,x);
y2=conv(x,h);
ny=[0:length(y1)-1];

subplot(2,2,1); stem(nx,x, 'd', 'filled')
xlabel('\itn'); ylabel('\itx\rm(\itn\rm)')
subplot(2,2,2); stem(nh,h, 'filled')
xlabel('\itn'); ylabel('\ith\rm(\itn\rm)')
subplot(2,2,3); stem(ny,y1, 'filled')
xlabel('\itn'); ylabel('\ity 1\rm(\itn\rm)')
subplot(2,2,4); stem(ny,y2, 'filled')
xlabel('\itn'); ylabel('\ity 2\rm(\itn\rm)')