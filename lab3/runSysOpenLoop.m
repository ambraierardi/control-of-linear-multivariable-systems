function [u_i, xtot] = runSysOpenLoop(A,B,T,x0,xf)

    syms tau real

    % define matrices
    eA_T = expm(A*T);
    eA_T_tau = expm(A*(T-tau));
    dG_cont = eA_T_tau*B*B'*eA_T_tau';
    G_cont = int(dG_cont,tau,0,T);

    iG_cont = inv(G_cont);

%x0 = [10 1]';
    u = B'*eA_T_tau'*iG_cont*(xf - eA_T*x0);
    u = simplify(u);
    
    dt = 0.1;
    time_F = [0:dt:T];
    time_F = time_F';
    LL = length(time_F);

    for j=1:length(u)
        uu = u(j);
        uu = symfun(uu,tau);
        u_i(j,:) = eval(uu(time_F));
    end


    xtot = zeros(length(x0),LL); % this is only needed to log the results and plot them
    xold = x0;

    for index = 1:LL-1

        dotx = A*xold + B*u_i(:,index);

        xnew = xold + dotx*dt;
        xtot(:,index+1) = xnew;

        xold = xnew;
    end

end