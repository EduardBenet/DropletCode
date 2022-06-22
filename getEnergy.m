function E = getEnergy(pos,sol,ContactAngle,param,A0)

a = param(1);
b = param(2);
n = param(5);
Lx = param(6);
p  = param(9); p = p/b;

y2 = pos + 1/sol(2) + 1/sol(2)*cos(sol(4));

A1 = 2*pi/sol(1)^2*(1+cos(sol(3)));
A2 = 2*pi/sol(2)^2*(1+cos(sol(4)));

Ac = 0; Vm = 0; Lt = 0;
for i = y2+(sol(5)-y2)/10000:(sol(5)-y2)/10000:sol(5)
    yu = i;
    yd = (i-(sol(5)-y2)/10000);
    xu = (1+p*i)*a*(1-abs(i/b)^n)^(1/n);
    xd = (1+p*(yd))*a*(1-abs(yd/b)^n)^(1/n);
    A = [xu,yu];
    B = [xd,yd];
    L = norm(A-B); Lt = Lt + L;
    Ac = Ac + L*2*pi*(Lx-(xu+xd)/2);
    Vm = Vm + pi*(Lx-(xu+xd)/2)^2*abs(i-(i-(sol(5)-y2)/1000));
end

E = A1+A2-Ac*cos(ContactAngle(b))-A0;