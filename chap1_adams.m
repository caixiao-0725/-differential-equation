function chap1_adams
% test the Adams interpolation method and the chap1_adams extrapolation method for ODE
% / u' = f(t,u),
% / u(0) = u0.

% by cheng xiao

u0=1;
T=2;
h = 0.1;
N = T/h;
t = 0:h:T
solu = exact1(t);


f = @f1;
u_inter_2s = adams_inter_2step(f,u0,t,h,N);
u_extra_2s = adams_extra_2step(f,u0,t,h,N);
figure(1)
plot(t,u_inter_2s,'*',t,u_extra_2s,'o',t,solu,'r')
legend('interpolation-2s','extrapolation-2s','exact')

u_inter_3s = adams_inter_3step(f,u0,t,h,N);
u_extra_3s = adams_extra_3step(f,u0,t,h,N);
figure(2)
plot(t,u_inter_3s,'*',t,u_extra_3s,'o',t,solu,'r')
legend('interpolation-3s','extrapolation-3s','exact')

end


function u = adams_inter_2step(f,u0,t,h,N)
% Adams interpolation method for ODE
u=zeros(1,N+1);
u(1)=u0;
%u(2)=u(1)+h*f(t(1),u(1));
u(2) = exact1(1*h);
eps_in = 1e-6;
K_in = 6;
for n = 2:N
    f_nm1 = f(t(n-1),u(n-1));
    f_n = f(t(n),u(n));
    s1 = u(n);
    du=1;
    k=1;
    while abs(du)>eps_in && k<=K_in
        s2 = u(n) + h/12*(5*f(t(n+1),s1)+8*f_n-f_nm1);
        du = s2-s1;
        s1 = s2;
        k=k+1;
    end
    u(n+1) = s2;
end

end

function u = adams_inter_3step(f,u0,t,h,N)
    % Adams interpolation method for ODE
    u=zeros(1,N+1);
    u(1)=u0;
    %u(2)=u(1)+h*f(t(1),u(1));
    %u(3)=u(2)+h*f(t(2),u(2));
    u(2) = exact1(1*h);
    u(3) = exact1(2*h);
    eps_in = 1e-6;
    K_in = 6;
    for n = 3:N
        f_nm1 = f(t(n-1),u(n-1));
        f_nm2 = f(t(n-1),u(n-2));
        f_n = f(t(n),u(n));
        s1 = u(n);
        du=1;
        k=1;
        while abs(du)>eps_in && k<=K_in
            s2 = u(n) + h/24*(9*f(t(n+1),s1)+19*f_n-5*f_nm1+f_nm2);
            du = s2-s1;
            s1 = s2;
            k=k+1;
        end
        u(n+1) = s2;
    end
    
end

function u = adams_extra_2step(f,u0,t,h,N)
% Adams extrapolation method for ODE
u=zeros(1,N+1);
u(1)=u0;
%u(2)=u(1)+h*f(t(1),u(1));
u(2) = exact1(1*h);
for n = 2:N
    u(n+1) = u(n) + h/2*(3*f(t(n),u(n))-f(t(n-1),u(n-1)));
end
end


function u = adams_extra_3step(f,u0,t,h,N)
% Adams extrapolation method for ODE
u=zeros(1,N+1);
u(1)=u0;
%u(2)=u(1)+h*f(t(1),u(1));
%u(3)=u(2)+h*f(t(2),u(2));
u(2) = exact1(1*h);
u(3) = exact1(2*h);
for n = 3:N
    u(n+1) = u(n) + h/12*(23*f(t(n),u(n))-16*f(t(n-1),u(n-1))+5*f(t(n-2),u(n-2)));
end
end

function f = f1(t,u)
f = -5*u;
end

function f = exact1(t)
f = exp(-5*t);
end


