%% Simple Electrical Circuit
% RLC circuit

clc
clear
close all

% define system parameters
R = 100;
Ca = 1e-3;
L = 1;

% define matrices
A = [0 1/Ca; -1/L -R/L];
B = [0 ;1/L];

C = [0 R];
D = 0;

% define time param
T = 1/100;
time = 2; 
N = round(time/T);

% define state vector, input, output
u = zeros(N,1);
vc = zeros(N,1);
i = zeros(N,1);
y = zeros(N,1);

u(1:N+1) = 5;
vc(1) = 0.0;
i(1) = 0.0;
xold = [vc(1) ; i(1)];
y(1) = 0.0;
yold = y(1);

t(1) = 0;

% integrate linear system
for index = 1:N
    xnew = xold + (A*xold + B*u(index))*T;
    ynew = yold + (C*xnew + D*u(index))*T;

    vc(index+1) = xnew(1);
    i(index+1) = xnew(2);
    y(index+1) = ynew;
    
    t(index+1) = t(index) + T;

    xold = xnew;
    yold = ynew;
end

%plot
figure;
subplot(2,1,1);
plot(t,vc,t,u,t,y); grid on;
title('Voltage')
ylabel('(volt)')
xlabel('Time(s)')
xlim([0,time])
legend('Vc','u','Vr')
hold on;

subplot(2,1,2);
plot(t,i,'g'); grid on;
title('Current')
ylabel('Current (A)')
xlabel('Time(s)')
xlim([0,time])

%svd is a command that tells you how far your matrix is from a singular
%matrix, even tough in an analitical way it would be a binary answer:
%either it is singular or not. 