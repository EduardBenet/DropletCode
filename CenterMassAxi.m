function z = CenterMassAxi(sol,pos,param)

a = param(1);
b = param(2);

n = param(5);
Lx = param(6);
p  = param(9); p = p/b;


R1 = 1/sol(1);
R2 = 1/sol(2);
alpha = sol(3);
beta = sol(4);
y2 = pos + 1/sol(2) + 1/sol(2)*cos(sol(4));

%Spherical Caps.
V1 = pi/(3*sol(1)^3)*(2+3*cos(sol(3))-cos(sol(3))^3);
V2 = pi/(3*sol(2)^3)*(2+3*cos(sol(4))-cos(sol(4))^3);

h1 = R1 + R1*cos(alpha);
h2 = R2 + R2*cos(beta);

z1 = 3*(2*R1-h1)^2/(4*(3*R1-h1))+sol(5) + 1/sol(1)*cos(sol(3));
z2 = -3*(2*R2-h2)^2/(4*(3*R2-h2))+pos+1/sol(2);

%Center
% Am = 0;xm = 0; 
ym = 0; Vm2 = 0;

% y2+(sol(5)-y2)/1000
% (sol(5)-y2)/1000
% sol(5)
for i = y2+(sol(5)-y2)/1000:(sol(5)-y2)/1000:sol(5)
    xu = (1+p*i)*a*(1-abs(i/b)^n)^(1/n);
    yu = i;
    yd = (i-(sol(5)-y2)/10000);
    xd = (1+p*(i-(sol(5)-y2)/1000))*a*(1-abs((i-(sol(5)-y2)/1000)/b)^n)^(1/n);
    ym  = ym + (yu+yd)/2*pi*(Lx-(xu+xd)/2)^2*(i-(i-(sol(5)-y2)/1000));
    Vm2 = Vm2 + pi*(Lx-(xu+xd)/2)^2*(i-(i-(sol(5)-y2)/1000));
end
% if Am < 0
%     Vm = 0;
% else
%     R = xm/Am;
%     Vm = 2*pi*R*Am;
% end

z = (V1*z1+V2*z2+ym)/(V1+Vm2+V2);