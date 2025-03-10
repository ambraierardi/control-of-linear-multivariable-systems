%% State Control and Stability Analysis

clc
clear
close all

% define system parameters
alfa = 30;
beta = 30;
gamma = 25;
m = 10;
Jw = 1.5;

% define matrices
A = zeros(6);
A(1:3,4:6) = diag(ones(1,3));
A(4:6,4:6) = diag([-alfa/m, -beta/m, -gamma/Jw]);
B = zeros(6,3);
B(4:6,1:3) = diag([1/m, 1/m, 1/Jw]);

C = zeros(3,6);
C(1:3,1:3) = diag([1,1,1]);

T = 30;
dt = 0.1;
N = T/dt;
time = 0:dt:T;
time_F = time;

%% Stabillity Analysis
% Set Kc (controlability)
Kc = zeros(3,6); % set the size of the matrix

% Set parameters for Kc
xi = 1;
wn = 3;

kc1 = -m*wn^2;
kc2 = -m*wn^2;
kc3 = -Jw*wn^2;
kc4 = alfa-2*m*xi*wn;
kc5 = beta-2*m*xi*wn;
kc6 = gamma-2*Jw*xi*wn;

Kc(1:3,1:3)=diag([kc1,kc2,kc3]);
Kc(1:3,4:6)=diag([kc4,kc5,kc6]);

% Set Ko (observability)
Ko = zeros(6,3); % set the size of the matrix

% Set parameters for Ko
xi = 1;
wn = 1;

kx1 = 2*xi*wn-alfa/m;
kx2 = 2*xi*wn-beta/m;
kx3 = 2*xi*wn-gamma/Jw;
kx4 = wn^2-kx1*alfa/m;
kx5 = wn^2-kx2*beta/m;
kx6 = wn^2-kx3*gamma/Jw;

Ko(1:3,1:3)=diag([kx1,kx2,kx3]);
Ko(4:6,1:3)=diag([kx4,kx5,kx6]);


%% Space state representation - Open Loop
% Initial condition and desired values

x0 = zeros(6,1);
xf = [5.0, 0.0, 0.0, 0.2, 0.0, 0.0]';

% Run Open Loop simulation
[u_d, xtot] = runSysOpenLoop(A,B,T,x0,xf);

% plot
% position
figure(1);
subplot(3,1,1)
title('Open Loop')
hold on;
plot(time_F, xtot(1,:)','k')
plot(time_F,ones(size(time_F))*xf(1),'r')
xlabel('Time [s]')
ylabel('Position [m]')
legend('x(t)','x_{ref}')
grid on;

% velocity
subplot(3,1,2);
plot(time,xtot(4,:)','g')
hold on;
plot(time,ones(size(time))*xf(4),'r')
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend('v(t)','v_{ref}')
grid on;

% acceleration
subplot(3,1,3);
plot(time, u_d(1,:)', time, u_d(2,:)',time, u_d(3,:)')
xlabel('Time [s]')
ylabel('Input [N]')
legend('F_x','F_y','N_w')
grid on;

xd = xtot; % save trajectory in xd vector

%% Closed Loop 
% without an observer (direct access on states)

syms tau
% start from a different initial state
x0 = [1.0, 3.0, 0.2, 3.5, 1.0, 0.1]';
xold = x0;
x_d = xtot;
for i=1:N+1
    % compute input u
    e=xd(:,i)-xold;
    u_cont(:,i) =u_d(:,i)-Kc*e;

    % compute the dynamic of the system dX
    dx = ((A*xold+B*u_cont(:,i))*dt);
    
    % compute the state X
    xnew = xold+dx;
    xtot_cont(:,i) = xnew; % save state vector
    % compute output Y
    Y(:,i) = C * xnew;
    
    xold = xnew;
end 


% plot 
% position x
figure(2);
subplot(3,2,1);
hold on;
plot(time, xtot_cont(1,:)','k')
plot(time, xd(1,:)','r')
xlabel('Time [s]')
ylabel('Position [m]')
legend('rx(t)','rx_{ref}')
grid on;

% velocity x
subplot(3,2,2);
plot(time, xtot_cont(4,:)','g')
hold on;
plot(time, xd(4,:)','r')
xlabel('Time [s]')
ylabel('surge [m/s]')
legend('Vx(t)','Vx_{ref}')
grid on;

% position y
subplot(3,2,3);
hold on;
plot(time, xtot_cont(2,:)','k')
plot(time, xd(2,:)','r')
xlabel('Time [s]')
ylabel('Position [m]')
legend('ry(t)','ry_{ref}')
grid on;

% velocity y
subplot(3,2,4);
plot(time, xtot_cont(5,:)','g')
hold on;
plot(time, xd(5,:)','r')
xlabel('Time [s]')
ylabel('sway [m/s]')
legend('Vy(t)','Vy_{ref}')
grid on;

% angle theta
subplot(3,2,5);
hold on;
plot(time, xtot_cont(3,:)','k')
plot(time, xd(3,:)','r')
xlabel('Time [s]')
ylabel('Angle [rad]')
legend('theta(t)','theta_{ref}')
grid on;

% angular velocity w
subplot(3,2,6);
plot(time, xtot_cont(6,:)','g')
hold on;
plot(time, xd(6,:)','r')
xlabel('Time [s]')
ylabel('Ang Velocity [rad/s]')
legend('w(t)','w_{ref}')
grid on;

%% Observer
% check Observability
syms t tau real
eA_t = expm(A*t);

% compute Gramiana
dG_obs = eA_t'*C'*C*eA_t;
G_obs = int(dG_obs,t,0,tau);
rank(G_obs);
disp(['rank(G_obs) = ' num2str(rank(G_obs))]);
if(rank(G_obs)==length(G_obs))
    disp('G_obs is full rank');
else
    disp('G_obs is NOT full rank');
end

%% Simulation 

% set standard deviation of the noise
r = 3;

x0 = [1.0, 3.0, 0.2, 0.5, 1.0, 0.1]'; % different initial states
xold = x0;

x0_est = zeros(6,1);

% estimated state and output
xold_est = x0_est;
Yest(min(size(Y)),1) = zeros(min(size(Y),1));

for i=1:N+1
    % compute input u
    %u_cont(:,i) = 

    % compute the dynamic of the system dx
    %dx = 

    % compute the state x
    %xnew = 
    xtot_cont(:,i) = xnew; % save the state vector for plotting

    % compute thhe output y
    %Y(:,i) = 
    xold = xnew;
    
    % generate a disturbed output 
    %Y_noise(:,i) =

    % compute input u for observer equation
    u_obs(:,i) = u_cont(:,i);

    % compute the dynamic of the Estimated system dx_est
    %dx_est = 

    % compute the Estimated state x_est
    %xnew_est = 
    xtot_obs(:,i) = xnew_est;  % save the state vector for plotting

    % compute the Estimated output y_est
    %Yest(:,i+1) = 
    xold_est = xnew_est;
end

% plot 
% position x
figure(3);
title('Closed Loop - Observer')
hold on;

subplot(3,2,1);
hold on;
plot(time_F, xtot_obs(1,:)','k')
plot(time_F, xd(1,:)','r')
xlabel('Time [s]')
ylabel('Position [m]')
legend('rx(t)','rx_{ref}')
grid on;

% velocity x
subplot(3,2,2);
plot(time_F, xtot_obs(4,:)','g')
hold on;
plot(time_F, xd(4,:)','r')
xlabel('Time [s]')
ylabel('surge [m/s]')
legend('Vx(t)','Vx_{ref}')
grid on;

% position y
subplot(3,2,3);
hold on;
plot(time_F, xtot_obs(2,:)','k')
plot(time_F, xd(2,:)','r')
xlabel('Time [s]')
ylabel('Position [m]')
legend('ry(t)','ry_{ref}')
grid on;

% velocity y
subplot(3,2,4);
plot(time_F, xtot_obs(5,:)','g')
hold on;
plot(time_F, xd(5,:)','r')
xlabel('Time [s]')
ylabel('sway [m/s]')
legend('Vy(t)','Vy_{ref}')
grid on;

% angle theta
subplot(3,2,5);
hold on;
plot(time_F, xtot_obs(3,:)','k')
plot(time_F, xd(3,:)','r')
xlabel('Time [s]')
ylabel('Angle [rad]')
legend('theta(t)','theta_{ref}')
grid on;

% angular velocity w
subplot(3,2,6);
plot(time_F, xtot_obs(6,:)','g')
hold on;
plot(time_F, xd(6,:)','r')
xlabel('Time [s]')
ylabel('Ang Velocity [rad/s]')
legend('w(t)','w_{ref}')
grid on;