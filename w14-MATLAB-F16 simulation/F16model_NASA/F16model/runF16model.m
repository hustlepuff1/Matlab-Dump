%=====================================================
%     Matlab Script File used to run the 
%     non-linear F-16 Simulation. The results
%     will also be plotted.
%
%   Author: Richard S. Russell
%
%  File made suitable for F-16 S-function dynamics.
%  - The simulink file "F16_trim.mdl" is used together
%    with the m-files "trim_fun.m" and "tgear.m"
%    to trim the aircraft.
%  - The simulink file "F16_openloop" is used to 
%    simulate the open loop response of the F-16 model.
%
%  L. Sonneveldt May 2006
%
%=====================================================

close all;
clear all;
clc;

%create the binary
mex F16_dyn.c

global altitude velocity fi_flag_Simulink

%% Ask user which simulation to run.
%%
newline = sprintf('\n');
disp('Which model would you like to use to trim the aircraft:')
disp('  1. Low Fidelity F-16 Trim')
disp('  2. High Fidelity F-16 Trim')
fi_flag = input('Your Selection:  ');
disp(newline);
disp(newline);

%% Determine from flag the correct simulation.
%%
if fi_flag == 1;
  fi_type = 'lofi';
  fi_flag_Simulink = 0;
elseif fi_flag == 2;
  fi_type = 'hifi';
  fi_flag_Simulink = 1;
else
  disp('Invalid selection');
  % break;
end

%% Trim aircraft to desired altitude and velocity
%%
altitude = input('Enter the altitude for the simulation (m)  :  ');
velocity = input('Enter the velocity for the simulation (m/s):  ');

%% Initial Conditions for trim routine.
%% The following values seem to trim to most
%% flight condition.  If the F16 does not trim
%% Change these values.
beta = 0;         % -
elevator = 0*pi/180;       % elevator, rad
alpha = 10*pi/180;           % AOA, rad
rudder = 0;         % rudder angle, rad
aileron = 0;         % aileron, rad
dth = 0.2;

% Initial Guess for free parameters
UX0 = [beta; elevator; alpha; aileron; rudder; dth];

% Initializing optimization options and running optimization:
OPTIONS = optimset('TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',5e+04,'MaxIter',1e+04);

iter = 1;
while iter == 1
    
    feval('F16_trim', [], [], [], 'lincompile');
    load_system('F16_trim');
   
    [UX,FVAL,EXITFLAG,OUTPUT] = fminsearch('trim_fun',UX0,OPTIONS);
   
    [cost, Xdot, xu, uu] = trim_fun(UX);
    
    disp('Trim Values and Cost:');
    disp(['cost   = ' num2str(cost)])
    disp(['dth    = ' num2str(uu(1)) ' -'])  
    disp(['elev   = ' num2str(uu(2)*180/pi) ' deg'])
    disp(['ail    = ' num2str(uu(3)*180/pi) ' deg'])
    disp(['rud    = ' num2str(uu(4)*180/pi) ' deg'])
    disp(['alpha  = ' num2str(xu(3)*180/pi) ' deg'])
    disp(['dLEF   = ' num2str(uu(5)*180/pi) ' deg'])
    disp(['Vel.   = ' num2str(xu(1)) ' m/s']) 
    disp(['pow    = ' num2str(xu(14)) ' %']) 
    flag = input('Continue trim rountine iterations? (y/n):  ','s'); 
    if flag == 'n'
        iter = 0;
    end
    feval('F16_trim', [], [], [], 'term'); clear F16_dyn;
    UX0 = UX;
end

% For simulink:
init_x = xu(1:14);
init_u = uu(1:4);
init_dlef = uu(5);

%sim open loop F-16 model
[t,x,y] = sim('F16_openloop',[0 10],simset('Solver','ode3','FixedStep',0.02)); clear F16_dyn;

%plots
figure(1);
subplot(3,1,1);plot(t,y(:,1));ylabel('Vt (m/s)' );
subplot(3,1,2);plot(t,y(:,2)*180/pi);ylabel('beta (deg)');
subplot(3,1,3);plot(t,y(:,3)*180/pi);ylabel('alpha (deg)');xlabel('t(s)');

figure(2);
subplot(4,1,1);plot(t,y(:,4));ylabel('q0' );
subplot(4,1,2);plot(t,y(:,5));ylabel('q1');
subplot(4,1,3);plot(t,y(:,6));ylabel('q2');
subplot(4,1,4);plot(t,y(:,7));ylabel('q3');xlabel('t(s)');

figure(3);
subplot(3,1,1);plot(t,y(:,8)*180/pi);ylabel('p (deg/s)' );
subplot(3,1,2);plot(t,y(:,9)*180/pi);ylabel('q (deg/s)');
subplot(3,1,3);plot(t,y(:,10)*180/pi);ylabel('r (deg/s)');xlabel('t(s)');

figure(4);
subplot(4,1,1);plot(t,y(:,11));ylabel('x (m)' );
subplot(4,1,2);plot(t,y(:,12));ylabel('y (m)');
subplot(4,1,3);plot(t,y(:,13));ylabel('z (m)');
subplot(4,1,4);plot(t,y(:,14));ylabel('power (%)');xlabel('t(s)');
      

