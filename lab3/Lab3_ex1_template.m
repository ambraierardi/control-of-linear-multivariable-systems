%% State Control and Stability Analysis

clc
clear
close all

% define system parameters
alfa = 30;
beta = 30;
gamma = 25;
m = 10;
Jw = 0.5;

% define matrices
A = zeros(6);
B = zeros(6,3);
C = zeros(3,6);
A(1:3,4:6)=eye(3);
A(4,4)=-alfa/m;
A(5,5)=-beta/m;
A(6,6)=-gamma/Jw;
B(4,1)=1/m;
B(5,2)=B(4,1);
B(6,3)=1/Jw;
C(1:3,1:3)=eye(3);

T = 50;
dt = 0.1;
N = T/dt;
time = 0:dt:T;

%% check Observability

syms tau
% compute Gramian
 eA_T = expm(A*T);
 eA_T_tau = expm(A*(T-tau));
 dG_obs = eA_T_tau'*C'*C*eA_T_tau;
 G_obs = int(dG_obs,tau,0,T);
 % check if Gramian matrice is full rank
 disp(['rank(G) = ' num2str(rank(G_obs))])
 if rank(G_obs)== min(size(G_obs))
    iG_cont = inv(G_obs);
 end
 
% print rank and check observability
rank(G_obs);

%% Stabillity Analysis
% find Kc and Ko
tau=3;
k1 = -alfa/m + 1/tau;
k2 = -beta/m + 1/tau;
k3 = -gamma/Jw + 1/tau;
K = [0 0 0 k1 0 0;
     0 0 0 0 k2 0;
     0 0 0 0 0 k3];

%% Space state representation
% set Initial condition and desired values

%X_0 = 
%X_d = X_desired
%u_d = u_desired

%X_old = 

% Loop for aqcuire Y
for i=0:N-1
    %u = 
    %dX = 
    %X_new = 
    %Y(i,:) = 
    X_old = X_new;
end
%% Observer

%Xest_0 = 
%Xest_old = 

Yest(0,:) = zeros(3,N);
for i=0:N-2
    %u = 
    %dXest = 
    %Xest_new = 
    %Yest(i+1,:) = 
    Xest_old = Xest_new;
    
end

%% Plot the result
% plot the state with respect to time (t,x)
% plot the cartesian position (rx,ry)













