clc;
clear;
digit=input('Enter 4 Digit\n','s');


p=(1/100)*str2double(digit(1));
m=(1/10)*str2double(digit(2));
xx=(1/100)*str2double(digit(3:4));

a0=0.2969; a1=-0.1260; a2=-0.3516; a3=0.2843; a4=-0.1015;

nx=60;
x=(1-cos(linspace(0,pi,nx)))./2;

Camber=[];

for i=1:size(x,2)
    if x(i)<m
        Camber(i)=(p/m^2)*((2*m*x(i))-(x(i))^2);
    elseif x(i)==m
        Camber(i)=p;
    elseif x(i)>m
        Camber(i)=(p/(1-m)^2)*((1-2*m)+2*m*x(i)-(x(i))^2);
    end
end

Thickness=(xx/0.2).*(a0*(x).^(1/2)+a1*(x)+a2*(x).^2+a3*(x).^3+a4*(x).^4); %Thickness

Delta=[0];
for i=1:size(x,2)
    if x(i)<m
        Delta(i)=atan(2*p/m^2)*(m-x(i));
    elseif x(i)>m
        Delta(i)=atan(2*p/(1-m)^2)*(m-x(i));
    end
end
    
xu=(x-(Thickness.*sin(Delta)));
yu=(Camber+(Thickness.*cos(Delta)));
xl=x+(Thickness.*sin(Delta));
yl=(Camber-(Thickness.*cos(Delta)));

nl=10;
l=linspace(0,2,nl);

sxu=zeros(nx,nl);
syu=zeros(nx,nl);
sxl=zeros(nx,nl);
syl=zeros(nx,nl);
scam=zeros(nx,nl);
sz=zeros(nx,nl);

for i=1:nl
    sxu(:,i)= xu;
    syu(:,i)=yu;
    sxl(:,i)= xl;
    syl(:,i)= yl;
    scamx(:,i)=x;
    scam(:,i)=Camber;
    sz(:,i)=l(i);
end

hold on
p1=surf(sxu,sz,syu);
p2=surf(sxl,sz,syl);

p3=surf(scamx,sz,scam);
legend('upper surface','lower surface','camberline')
p1.FaceColor='r';
p2.FaceColor='y';
p3.FaceColor='b';
axis([0 1 0 2 -0.05 0.20]);
alpha(0.7);
grid on
view([30 30]);
