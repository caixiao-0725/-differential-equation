function chap1_euler
% Euler method for solving ODEs
% u'=f(t,u), u(0)=u0
% f(t,u) = -5u

% by cheng xiao
T=1;
h=0.01;
t=0:h:T;
N=length(t)-1;

%精确解
solu = exp(-5*t);

u0 =1;
f = @f1;
u_euler = euler(f,u0,t,h,N);
u_implicit_euler = implicit_euler(f,u0,t,h,N);
u_modified_euler = modified_euler(f,u0,t,h,N);

figure(1)
plot(t,solu,'r',t,u_euler,'*b',t,u_implicit_euler,'gs',t,u_modified_euler,'s')
legend('exact','euler','implicit euler','modified euler')
end


function u = euler(f,u0,t,h,N)
u = zeros(N+1,1);
u(1) = u0;
for i=1:N
    u(i+1) = u(i) + h*f(t(i),u(i));
end
end

function u = implicit_euler(f,u0,t,h,N)
u = zeros(N+1,1);
u(1) = u0;
eps_in = 1e-6;
K_in = 6;
for n=1:N
    s1 = u(n);
    k=1;
    du =1;
    % du 就是看两个迭代值之间的差值是不是足够小
    while k<=K_in && abs(du)>eps_in
        s2 = u(n)+h*f(t(n+1),s1);
        du = s2-s1;
        s1 = s2;
        k = k+1;
    end
    u(n+1) = s1;
end
end

function u = modified_euler(f,u0,t,h,N)
u = zeros(N+1,1);
u(1) = u0;
eps_in = 1e-6;
K_in = 6;
for n=1:N
    fn = f(t(n),u(n));
    s1 = u(n);
    k=1;
    du = 1;
    while k<=K_in && abs(du)>eps_in
        s2 = u(n)+h/2*(fn+f(t(n+1),s1));
        du = s2-s1;
        s1 = s2;
        k = k+1;
    end
    u(n+1) = s1;
end
end

function f = f1(t,u)
f=-5*u;
end
