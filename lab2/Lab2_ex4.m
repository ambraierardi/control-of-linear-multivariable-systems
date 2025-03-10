%% Adding noise on the output y
% changing standard deviation

r = [0.0 1.0 3.0 9 27 100]; % standard deviation array
clear y_noise y_f_val x0_est_tot y_f

for j=1:length(r)

    y_noise(1,:) = xtot(1,:) + r(j).*randn(1,length(time_F));

    dy_f = eA_t_tau*B*u_fun;
    y_f = C*int(dy_f,tau,0,t);

    y_f_fun1 = symfun(y_f(1),t);
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

    x0_est = iG_obs_val*znew;
    x0_est_tot(:,j) = x0_est;

end

% plot
figure(3);
%subplot(2,1,1);
plot(r,x0_est_tot(1,:),'o')
hold on
grid minor
plot(r,ones(size(r))*x0(1))
xlabel('Output standard deviation [m]')
ylabel('distance [m]')
legend('r_{est}','r_{true}')



