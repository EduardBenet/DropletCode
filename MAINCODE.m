close all; clc;
options = optimset('TolFun',1e-16,'TolX',1e-16,'MaxFunEvals',1200,'display','off');

%Droplet going through a pore.
%Code by E.Benet in F.Vernerey group.

%Input parameters are the real values, the code then non-dimensinoalizes
%them, solves the problem, and returns a non-dimensional solution

%The problem starts (and ends) with a vesicle outside the pore, which means that it
%might take a while until it founds a plausible solution, and it might give
%non-physical solutions ones the vesicle gets out.

%Since the problem works controlling the bottom edge of the vesicle the
%second part of the results will be less accurate than the first part.

%% Input Parameters
a_real  = 2;   %Geometry Width
b_real  = 2;   %Geometry Height
n       = 2;   %Geometry edge sharpness, large values (above 100) might not give a good solution.
p       = 0;   %Geometry Slope (0 to 1)
d_real  = 1;   %Pore smallest diameter (keep at 1 do not change)
R0_real = 1;   %Vesicle Radius (1 means the vesicle is twice the diameter)
gamma   = 0.015; %Surface Tension
ContactAngle = @(y) pi;                %Contact Angle constant
plott = 1; %Plot (1) or not(0)


%% Non-Dimensionalize problem
a  = a_real/d_real;
b  = b_real/d_real;
R0 = R0_real/d_real;
d  = 1;
theta = ContactAngle(b);

[x0,y0] = InitialPoint(n,p,a,b,ContactAngle);
if (n == 2) && (a == b)
    Ip    = sqrt((a+R0)^2-(a+d/2)^2);
    if theta == pi
        pos = sqrt(-(a+d/2)^2+(a+R0)^2)-R0;
    else
        pos = a*sin(theta-pi/2);
    end
    Range = pos:-0.01/d_real:-b-2*R0;
else
    END = R0-R0*cos(ContactAngle(-b));
    pos   = y0; %Initial position of the bottom edge of the vesicle, careful because this value is not arbitrary and the code might fail
    Range = pos:-0.01/d_real:-b-2*R0;
    Ip = b;
end

A0 = 4*pi*R0^2;
V0 = 4/3*pi*R0^3;

%% Code
%Find Edge of the geometry (maximum x)
ymax = GoldenSearch(a,b,n,p);
xmax = (1+p*ymax/b).*(1-abs(ymax/b).^n).^(1/n)*a;

Lx = xmax+d/2;
param = [a,b,Range(2),gamma,n,Lx,R0,d,p];

%Define and ploy geometry
y = -b:0.01:b;
x = (1+p*y/b).*(1-abs(y/b).^n).^(1/n)*a;

if plott == 1
    subplot(2,2,1);
    fig = figure(1); hold on; axis equal;
    set(fig,'units','normalized','outerposition',[0 0 1 1]);
    
    angle = linspace(0,2*pi);
    X = (x-Lx)'*cos(angle);
    Y = (x-Lx)'*sin(angle);
    Z = y'*ones(size(angle));
    surf((X+Lx)*d_real,Y*d_real,Z*d_real,zeros(size(Z)),'LineWidth',2,'FaceAlpha',.5);
    shading interp;
    light;
    lighting gouraud;
    camlight left;
    view(10,20)
    axis([0 (2*a*(1+p)+d)*d_real -d_real d_real (-R0-b)*d_real (R0+b)*d_real])
end
y = -b:0.0001:b;
x = (1+p*y/b).*(1-abs(y/b).^n).^(1/n)*a;

%Initial Guess;
thetaC = pi;
U0 = [1/(R0); 2/R0; pi/6+(pi-thetaC); 0.8*pi+(pi-thetaC);  b; 0];

UB = [2/d 2/d 2*pi 2*pi b inf];
LB = [0, 0,0,0,-b,-inf];

sol = lsqnonlin(@(U) Residual3D(U,param,ContactAngle),U0,LB,UB,options);

if isreal(sol) == 0
    if isreal(sol(1)) == 0
        sol(1) = 2/R0;
    end
    if isreal(sol(2)) == 0
        sol(2) = 1/R0;
    end
    if isreal(sol(3)) == 0
        sol(3) = pi/2;
    end
    if isreal(sol(4)) == 0
        sol(4) = pi/2;
    end
    if isreal(sol(5)) == 0
        sol(5) = b;
    end
    if isreal(sol(6)) == 0
        sol(6) = 0;
    end
end

%Solve using the previous solution as first guess
Position = 0;CritPres = 0; E = 0;
h1 = []; h2 = []; h3 = []; h4 = []; h5 = []; h6 = []; h7 = []; h8 = [];
count = 0;
y2_old = 0.95*b;
disp('Refine End')
Inc0 = -b/100; Inc = -0.1;Bcount = 100;Icount = 0;
while abs(Inc0)>1e-10
    if y2_old<-0.95*b && abs(Inc)>b/1e4
        Inc = -b/1e4;
    end
    y2 = y2_old + Inc;
    if y2 < -b
        Inc0 = Inc0/10;Inc = Inc0;
        Bcount = 0;
    else
        param(3) = y2;
        LB = [0, 0,0,0,-b,-inf];
        sol = lsqnonlin(@(U) Residual3D2(U,param,ContactAngle),sol,LB,UB,options);
        
        pos = y2 - 1/sol(2) - 1/sol(2)*cos(sol(4));
        y1 = sol(5);
        tol = norm(Residual3D2(sol,param,ContactAngle));
        param(3) = pos;
        
        if abs(y2) > b || (y2 - y1)>1e-5 || tol > 1e-5
            Bcount = Bcount + 1;
            disp(['Non Physical Solution. N = ',num2str(n),' p = ',num2str(p),' R = ',num2str(R0_real),' Inc: ',num2str(Inc),'  y2:',num2str(y2_old)])
            if Bcount < 100
                Inc = Inc/10;
            else
                y2_old = y2;
            end
            sol = U0;
            if Bcount > 5 && Bcount <100
                break
            end
        else
            U0 = sol;
            count = count + 1;Bcount = 0;Icount = Icount + 1;
            if count == 1
                first_y = y2;
                first_sol = U0;
                Inc0 = Inc0/10; Inc = Inc/10;
            end
            Position(count) = CenterMassAxi(sol,pos,param)*d_real/b;
            CritPres(count) = sol(6);
            E(count) = getEnergy(pos,sol,ContactAngle,param,A0);
            
            if plott == 1
                %Plot 3D
                delete(h1); delete(h2); delete(h3); delete(h4); delete(h5); delete(h6); delete(h7); delete(h8);
                x2 = (1+p*y2/b)*(1-abs(y2/b)^n)^(1/n)*a;
                x1 = (1+p*y1/b)*(1-abs(y1/b)^n)^(1/n)*a;
                
                Center1 = [xmax+d/2,sol(5) + 1/sol(1)*cos(sol(3))];
                Center2 = [xmax+d/2,pos+1/sol(2)];
                pplot1 = 1/sol(1)*cos((pi/2):0.001:(3*pi/2-sol(3)));
                pplot2 = 1/sol(1)*sin((pi/2):0.001:(3*pi/2-sol(3)));
                pplot3 = 1/sol(2)*cos((pi/2+sol(4)):0.001:(3*pi/2));
                pplot4 = 1/sol(2)*sin((pi/2+sol(4)):0.001:(3*pi/2));
                
                pplot6 = linspace(y1,y2,100);
                pplot5 = (1+p*pplot6/b).*(1-abs(pplot6/b).^n).^(1/n)*a-Lx;
                thetaC1 = ContactAngle(y1);
                thetaC2 = ContactAngle(y2);
                
                subplot(2,2,1);
                X = [pplot1'*cos(angle)+Center1(1);pplot5'*cos(angle)+Center1(1);pplot3'*cos(angle)+Center2(1)]; Y = [pplot1'*sin(angle);pplot5'*sin(angle);pplot3'*sin(angle)]; Z = [pplot2'*ones(size(angle))+Center1(2);pplot6'*ones(size(angle));pplot4'*ones(size(angle))+Center2(2)];
                h1 = surf(X*d_real,Y*d_real,Z*d_real,ones(size(Z)),'LineWidth',2,'FaceAlpha',1);
                shading interp;
                axis([0 (2*a*(1+p/b)+d) -1 1 -R0-b R0+b]*d_real)
                
                %Plot Profile
                subplot(2,2,3); hold on;
                h7 = plot([pplot1,pplot5,pplot3]*d_real,[pplot2+Center1(2),pplot6,pplot4+Center2(2)]*d_real,'r','LineWidth',2);
                h3 = text((x2-Lx)*d_real,y2*d_real,['   ',num2str(thetaC2*180/pi)]);
                h4 = text((x1-Lx)*d_real,y1*d_real,['   ',num2str(thetaC1*180/pi)]);
                h5 = plot((x1-Lx)*d_real,y1*d_real,'ro');
                h6 = plot((x2-Lx)*d_real,y2*d_real,'ro');
                plot(x-Lx,y,'k')
                axis equal
                title('Contact Angle at each Location')
                %Plot Pressure
                subplot(2,2,2); hold on;
                set(gca,'XDir','reverse')
                h2 = plot(Position,CritPres,'b','LineWidth',2); hold on;
                xlabel('Position Center of mass')
                ylabel('Pressure')
                title(['Pressure at each position. Diff = ',num2str(max(CritPres)-abs(min(CritPres)))])
                %Plot Energy
                subplot(2,2,4)
                h8 = plot(Position,E,'b','LineWidth',2);
                set(gca,'XDir','reverse')
                ylabel('Energy')
                pause(0.1)
            end
            sol_old = sol;
            y2_old = y2;
            pos_old = pos;
            if abs(Inc) < b/10000 && Icount > 20
                Icount = 0;
                Inc = Inc0;
            end
        end
    end
end

disp('Refine Top')
y2_old = first_y-b/100; sol = first_sol; U0 = first_sol;
done = 0; Inc0 = b/100; Inc = 0.0001;Bcount = 100;Icount = 0;
while abs(Inc0)>1e-10
    if y2_old>0.95*b && abs(Inc)>b/1e4
        Inc = b/1e4;
    end
    y2 = y2_old + Inc;
    if y2 > b
        Inc0 = Inc0/10;Inc = Inc0;
        Bcount = 0;
    else
        param(3) = y2;
        LB = [0, 0,0,0,-b,-inf];
        sol = lsqnonlin(@(U) Residual3D2(U,param,ContactAngle),sol,LB,UB,options);
        pos = y2 - 1/sol(2) - 1/sol(2)*cos(sol(4));
        y1 = sol(5);
        tol = norm(Residual3D2(sol,param,ContactAngle));
        param(3) = pos;
        if abs(y2) > b || (y2 - y1)>1e-5 || tol > 1e-10
            Bcount = Bcount + 1;
            disp(['Non Physical Solution ',num2str(Bcount),'  N = ',num2str(n),' p = ',num2str(p),' R = ',num2str(R0_real),' Inc: ',num2str(Inc),'  y2:',num2str(y2)])
            Inc = Inc/10;
            sol = U0;
            if Bcount > 5 && Bcount <100
                break
            end
        else
            U0 = sol;
            Bcount = 0;Icount = Icount + 1;
            
            Position = [CenterMassAxi(sol,pos,param)*d_real/b,Position];
            CritPres = [sol(6),CritPres];
            E = [getEnergy(pos,sol,ContactAngle,param,A0),E];
            
            if plott == 1
                %Plot 3D
                delete(h1); delete(h2); delete(h3); delete(h4); delete(h5); delete(h6); delete(h7); delete(h8);
                x2 = (1+p*y2/b)*(1-abs(y2/b)^n)^(1/n)*a;
                x1 = (1+p*y1/b)*(1-abs(y1/b)^n)^(1/n)*a;
                
                Center1 = [xmax+d/2,sol(5) + 1/sol(1)*cos(sol(3))];
                Center2 = [xmax+d/2,pos+1/sol(2)];
                pplot1 = 1/sol(1)*cos((pi/2):0.001:(3*pi/2-sol(3)));
                pplot2 = 1/sol(1)*sin((pi/2):0.001:(3*pi/2-sol(3)));
                pplot3 = 1/sol(2)*cos((pi/2+sol(4)):0.001:(3*pi/2));
                pplot4 = 1/sol(2)*sin((pi/2+sol(4)):0.001:(3*pi/2));
                
                pplot6 = linspace(y1,y2,100);
                pplot5 = (1+p*pplot6/b).*(1-abs(pplot6/b).^n).^(1/n)*a-Lx;
                thetaC1 = ContactAngle(y1);
                thetaC2 = ContactAngle(y2);
                
                subplot(2,2,1);
                X = [pplot1'*cos(angle)+Center1(1);pplot5'*cos(angle)+Center1(1);pplot3'*cos(angle)+Center2(1)]; Y = [pplot1'*sin(angle);pplot5'*sin(angle);pplot3'*sin(angle)]; Z = [pplot2'*ones(size(angle))+Center1(2);pplot6'*ones(size(angle));pplot4'*ones(size(angle))+Center2(2)];
                h1 = surf(X*d_real,Y*d_real,Z*d_real,ones(size(Z)),'LineWidth',2,'FaceAlpha',1);
                shading interp;
                axis([0 (2*a*(1+p/b)+d) -1 1 -R0-b R0+b]*d_real)
                
                %Plot Profile
                subplot(2,2,3); hold on;
                h7 = plot([pplot1,pplot5,pplot3]*d_real,[pplot2+Center1(2),pplot6,pplot4+Center2(2)]*d_real,'r','LineWidth',2);
                h3 = text((x2-Lx)*d_real,y2*d_real,['   ',num2str(thetaC2*180/pi)]);
                h4 = text((x1-Lx)*d_real,y1*d_real,['   ',num2str(thetaC1*180/pi)]);
                h5 = plot((x1-Lx)*d_real,y1*d_real,'ro');
                h6 = plot((x2-Lx)*d_real,y2*d_real,'ro');
                plot(x-Lx,y,'k')
                axis equal
                title('Contact Angle at each Location')
                %Plot Pressure
                subplot(2,2,2); hold on;
                set(gca,'XDir','reverse')
                h2 = plot(Position,CritPres,'b','LineWidth',2); hold on;
                xlabel('Position Center of mass')
                ylabel('Pressure')
                title(['Pressure at each position. Diff = ',num2str(max(CritPres)-abs(min(CritPres)))])
                %Plot Energy
                subplot(2,2,4)
                h8 = plot(Position,E,'b','LineWidth',2);
                set(gca,'XDir','reverse')
                ylabel('Energy')
                pause(0.1)
            end
            sol_old = sol;pos_old = pos;
            y2_old = y2;
            if abs(Inc) < b/10000 && Icount > 20
                Icount = 0;
                Inc = Inc0;
            end
        end
    end
    if Inc <1e-20
        break
    end
end
