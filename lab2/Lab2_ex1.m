%% Study the Controllability

clc
clear
close all

% define system parameters
k = 1e-02;
% define matrices
A = [0 1; 0 -k];
B = [0 1]';
T = 50;
dt = 0.1;
xf = [-20 +2]'; % reference states

syms tau

% define matrices
eA_T = expm(A*T);
eA_T_tau = expm(A*(T-tau));
dG_cont = eA_T_tau*B*B'*eA_T_tau';
G_cont = int(dG_cont,tau,0,T);

% check if Gramian matrice is full rank
disp(['rank(G) = ' num2str(rank(G_cont))])

if rank(G_cont)== min(size(G_cont))
    iG_cont = inv(G_cont);

    x0 = [10 1]'; 
    u = B'*eA_T_tau'*iG_cont*(xf - eA_T*x0);
    u = simplify(u);

    u = symfun(u,tau);

    time_F = [0:dt:T];
    time_F = time_F';
    LL = length(time_F);

    uvalue = eval(u(time_F));

    
    xtot = zeros(2,LL); % this is only needed to log the results and plot them
    xold = x0;
    for index = 1:LL-1

        dotx = A*xold + B*uvalue(index);

        xnew = xold + dotx*dt;
        xtot(:,index+1) = xnew;

        xold = xnew;
    end

else
    disp('G is not full rank');
end


% plot
% position
figure(1);
subplot(3,1,1);
plot(time_F,xtot(1,:)','k')
hold on;
plot(time_F,ones(size(time_F))*xf(1),'r')
xlabel('Time [s]')
ylabel('Position [m]')
legend('x(t)','x_{ref}')
grid on;

% velocity
subplot(3,1,2);
plot(time_F,xtot(2,:)','g')
hold on;
plot(time_F,ones(size(time_F))*xf(2),'r')
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend('v(t)','v_{ref}')
grid on;

% acceleration
subplot(3,1,3);
plot(time_F,uvalue','b')
xlabel('Time [s]')
ylabel('Input [m/s^2]')
grid on;







