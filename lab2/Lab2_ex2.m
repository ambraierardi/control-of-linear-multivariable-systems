%% Study the observability

%clc
%clear
%close all

% define system parameters
k = 1e-02;
% define matrices
A = [0 1; 0 -k]; 
B = [0 1]';

C = [1 0]; 
D = 0;

T = 50;
dt = 0.1;
time = 0:dt:T;
time = time';

syms t tau

% define matrices
eA_t_tau = expm(A*(t-tau));
eA_t = expm(A*t);

y_m(1,:) = xtot(1,:);

% curve fitting for input signal u
p = polyfit(time,uvalue,2);

% equation of input used for integration u(tau)
u_fun = [tau^2 tau 1]*p';

dy_f = eA_t_tau*B*u_fun;
y_f = C*int(dy_f,tau,0,t);

y_f_fun1 = symfun(y_f(1),t);
y_f_val(1,:) = eval(y_f_fun1(time));

eA_t_fun = symfun(eA_t,t);

y_l = y_m - y_f_val;
zold = [0, 0]';

for i=1:length(y_l)
    eA_t_val = eval(eA_t_fun(time(i)));
    dz = eA_t_val'*C'*y_l(:,i);
    znew = zold + dz*dt;
    zold = znew;
end

% compute Gramian
dG_obs = eA_t'*C'*C*eA_t;
G_obs = int(dG_obs,t,0,tau);
rank(G_obs)
if(rank(G_obs)==length(G_obs))
    iG_obs = inv(G_obs);
    iG_obs = simplify(iG_obs);
    iG_obs = symfun(iG_obs,tau);
    iG_obs_val = eval(iG_obs(T));

    x0_est = iG_obs_val*znew
else
    disp('G is not full rank');
end



