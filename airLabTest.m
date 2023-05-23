% Plot streamlines and pressure coefficent of an unviscious,uncomprimibile,irrotational
% flow around a cylinder section (radius = 1) that spins around the z axis (coming out
% of the xy plane).
% This result is achieved by superimposition of elementary solution of the potential
% function PHI, where [Ux , Uy] = GRAD(PHI) which comprehend Uniform Stream,Doublet,Vortex.
% This case rappresents a good example of Magnus Effect,that is the reason why spinning
% balls have an effected trajectory.And rappresent the basis of Kutta-Joukowski theory.
%
%
% INPUT
% V_i = Asymptotic Speed
% G = Angular Speed (positive anti-clockwise)
%
% -----EXAMPLE------
% V_i = 20
% G = 50
%
% Created by: Dario Isola
% DATA : 24/03/2004
% Politecnico di Milano - Aerospatial Engeeniering Departement
%----------
%Modification log
%Author Date Description
%----------------------------------------------------------------------------------
%Yogesh PARTE,IMT Toulouse 15 Jan 2010 Added Cp distribution plot
% over a cylinder see
% figure(4), removed for loop, added comments
%----------------------------------------------------------------------------------
close all;
clear all;
%% Input section
V_i = input(' Asymptotic speed V_0 [m/s] = ');
G = input(' Circulation Value G [rad/s] [Anti-clockwise] = ');
%% Actual computation
%Radius of the circle
a = 1 ;
c =-a*2;
b =a*2;
% Number of intervals
n =a*50;
[x,y]=meshgrid([c:(b-c)/n:b],[c:(b-c)/n:b]');
warning off;
%%Preliminary DATA & purification
% Set values of X and Y inside the cylinder to zero
[I J]=find( (x.^2+y.^2) < a);
if ~isempty(I)
x(I,J) = 0;
y(I,J) = 0;
end
%Definition of polar variables
rho=sqrt(x.^2+y.^2);
theta=atan2(y,x);
% Creation of the streamline function
z=V_i.*sin(theta).*rho.*(1-(a^2./(rho.^2)))-G*log(rho)/(2*pi);
%% Generate unite cicle for plotting
n=100;
r=ones(1,n+1)*a;
t=[0:2*pi/n:2*pi];
Xcircle = r.*cos(t);
Ycircle = r.*sin(t);
%% Plot the data
% Streamline Plot
figure(1);
contour(x,y,z,25);
colorbar;
hold on;
fill(Xcircle,Ycircle,'k');
title('Stream Lines');
xlabel('x \rightarrow');
ylabel('y \rightarrow');
axis square;
%% SECOND PART
% Reproduce streamlines
figure(2);
contour(x,y,z,15);
hold on;
fill(Xcircle,Ycircle,'k');
%Compute velocity at new X,Y coordinates for better viewing
%Creation of vectors around the circle
x=[-a*2:a/3:a*2];
[x]=meshgrid(x);
y=x';
% Set values of X and Y inside the cylinder to zero
[I J]=find( (x.^2+y.^2) < a);
if ~isempty(I)
x(I,J) = 0;
y(I,J) = 0;
end
r=sqrt(x.^2+y.^2);
theta=atan2(y,x);
ur=V_i*cos(theta).*(1-a^2./(r.^2));
ut=-V_i*sin(theta).*(1+a^2./(r.^2))+G./(2*pi*r);
u=ur.*cos(theta)-ut.*sin(theta);
v=ur.*sin(theta)+ut.*cos(theta);
%Vectors and Filled Circle plotting
hold on;
quiver(x,y,u,v);
title('Speed Vectors')
xlabel('x \rightarrow');
ylabel('y \rightarrow');
axis square;
grid off;
% Analytical solution for Cp distribution
t=0:.1:2*pi;
cp = 1 - 4*sin(t).^2 + 2* G / (pi*a*V_i) *sin(t) - (2* G/ (pi*a*V_i) )^2 ;
% Non lifting solution
cp_sim = 1 - 4*sin(t).^2 ;
% Lift from Kutta Jokowaski theorm
L = - 1.225*V_i*G;
L = strcat('Kutta Joukowski Lift: ',num2str(L),' [N]');
% plot cp distribution
figure(3);
plot(t,cp,t,cp_sim,'--r');
axis([0 2*pi min(cp) max(cp_sim)]);
title(strvcat('Pressure coefficient around the surface (standard air density)',L));
xlabel('Theta (angle with horizontal)')
ylabel('C_p')
legend('Lifting solution','Symmetrical solution');
grid on;
% Plot cp distribition on cylinder
figure(4);
%filled cirlcle
cpx = cp.*cos(t);
cpy = cp.*sin(t);
fill(Xcircle,Ycircle,'y');
hold on;
scale =2; % used to control length of Cp in quiver plot
quiver(a*cos(t),a*sin(t),cpx,cpy,scale);
title('Pressure coefficient around the cylinder surface (standard air density)');
xlabel('x \rightarrow');
ylabel('y \rightarrow');
ylim([-1.5 1.5]);
axis square;