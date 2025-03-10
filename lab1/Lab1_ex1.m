clc
clear
close all

% Simple Mechanical System
% mass-spring_friction system

% define system parameters
k = 100;
m = 1;
f = 5;

% define matrices
A = [0 1; -k/m -f/m];
B = [0; 1/m];

% define time param
T = 1/100;
time = 5; 
N = round(time/T);

% define state vector, input
u = zeros(N,1);
v = zeros(N,1);
r = zeros(N,1);
t = zeros(N,1);

% initial state and input
t(1) = 0;
u(1) = 2;
v(1) = 0;
r(1) = 0.2;
xold = [v(1) ; r(1)];

% integrate linear system
for index = 1:N
    xnew = xold + (A*xold + B*u(index))*T;
    v(index+1) = xnew(1);
    r(index+1) = xnew(2);
    t(index+1) = t(index) + T;

    xold = xnew;
end

% plot position and velocity
figure;
subplot(2,1,1);
plot(t,v,'g'); grid on;
ylabel('Velocity (m/s)')
xlabel('Time(s)')
hold on;

subplot(2,1,2);
plot(t,r,'r'); grid on;
ylabel('distance (m)')
xlabel('Time(s)')