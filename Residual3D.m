function Fk = Residual3D(U,param,ContactAngle)
a = param(1);
b = param(2);
P = param(3);
n = param(5);
Lx = param(6);
R0 = param(7);
p  = param(9); p = p/b;

if isreal(U(1)) == 0
    U(1) = 2/R0;
end
if isreal(U(2)) == 0
    U(2) = 1/R0;
end
if isreal(U(3)) == 0
    U(3) = pi/2;
end
if isreal(U(4)) == 0
    U(4) = pi/2;
end
if isreal(U(5)) == 0
    U(5) = b;
end
if isreal(U(6)) == 0
    U(6) = 0;
end

y2 = P + 1/U(2) + 1/U(2)*cos(U(4));
if y2 > b || y2 < -b || (isreal(y2)==0)
    y2 = max(P,-0.99*b);
end
x2 = (1+p*y2)*(1-abs(y2/b)^n)^(1/n)*a;
x1 = (1+p*U(5))*(1-abs(U(5)/b)^n)^(1/n)*a;

alphae = 3*pi/2-U(3);
betae  = p*(1-abs(U(5)/b).^n).^(1/n)*a+(1+p*U(5)).*(a/n*(1-abs(U(5)/b).^n).^(1/n-1)).*(-n*abs(U(5)/b).^(n-1)).*sign(U(5))/b;
alphai = U(4)+pi/2;
betai  = p*(1-abs(y2/b).^n).^(1/n)*a+(1+p*y2).*(a/n*(1-abs(y2/b).^n).^(1/n-1)).*(-n*abs(y2/b).^(n-1)).*sign(y2)/b;

Vm = 0;
for i = y2+(U(5)-y2)/1000:(U(5)-y2)/1000:U(5)
    xu = (1+p*i)*a*(1-abs(i/b)^n)^(1/n);
    xd = (1+p*(i-(U(5)-y2)/1000))*a*(1-abs((i-(U(5)-y2)/1000)/b)^n)^(1/n);
    Vm = Vm + pi*(Lx-(xu+xd)/2)^2*abs((U(5)-y2)/1000);
end

Ve = pi/(3*U(1)^3)*(2+3*cos(U(3))-cos(U(3))^3);
Vi = pi/(3*U(2)^3)*(2+3*cos(U(4))-cos(U(4))^3);

thetaC1 = ContactAngle(U(5));
thetaC2 = ContactAngle(y2);

Fk = zeros(6,1);

Fk(1) = U(6)-(U(2)-U(1));
Fk(2) = 4/3*pi*R0^3-Ve-Vm-Vi;
% atan2(b*sin(betae)^(2/n-1)*cos(betae),-a*cos(betae)^(2/n-1)*sin(betae))*180/pi
% atan2(-cos(alphae),sin(alphae))*180/pi
Fk(3) = pi-thetaC1-(atan2(1,betae)-atan2(-cos(alphae),sin(alphae)));
Fk(4) = pi-thetaC2+(atan2(1,betai)-atan2(-cos(alphai),sin(alphai)));
Fk(5) = sin(U(3))-(Lx-x1)*U(1);
Fk(6) = sin(U(4))-(Lx-x2)*U(2);
end



