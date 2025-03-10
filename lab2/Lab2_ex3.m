%% Adding noise on the output y
%close all
r = 2.0; % standard deviation

y_noise(1,:) = xtot(1,:) + r.*randn(1,length(time_F)); 
%randn gives a normally distributed error
%gaussian noise is not that bad, it's zero-mean6

dy_f = eA_t_tau*B*u_fun; 
y_f = C*int(dy_f,tau,0,t); %change to T

y_f_fun1 = symfun(y_f,t);
y_f_val(1,:) = eval(y_f_fun1(time_F));


eA_t_fun = symfun(eA_t,t);

y_l = y_noise - y_f_val;
zold = [0, 0]';
for i=1:max(size(y_l))
    eA_t_val = eval(eA_t_fun(time_F(i)));
    dz = eA_t_val'*C'*y_l(:,i);
    znew = zold + dz*dt;
    zold = znew;
end

x0_est = iG_obs_val*znew


% plot
% r distance
figure (2);
plot(time_F,y_noise(1,:)','k')
hold on;
xlabel('Time [s]')
ylabel('Position [m]')
grid on;
 
