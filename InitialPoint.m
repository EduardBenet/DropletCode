function [x,y] = InitialPoint(n,p,a,b,ContactAngle)

p = p/b;
% y = -b:0.01:b;y = y';
% x = (1+p*y).*(1-abs(y/b).^n).^(1/n)*a;
% figure
% plot(x,y,'LineWidth',2); hold on;
% tp = [p*(1-abs(y/b).^n).^(1/n)*a+(1+p*y).*(a/n*(1-abs(y/b).^n).^(1/n-1)).*(-n*abs(y/b).^(n-1)).*sign(y)/b,ones(size(y))];
% tp = [tp(:,1).*(tp(:,1).^2+tp(:,2).^2).^-0.5,tp(:,2).*(tp(:,1).^2+tp(:,2).^2).^-0.5];

tol = 1;iter = 0;
y_up = b;y_down = 0;
while abs(tol) >1e-12
    iter = iter + 1;
    y = (y_up+y_down)/2;
    angle = atan2(1,p*(1-abs(y/b).^n).^(1/n)*a+(1+p*y).*(a/n*(1-abs(y/b).^n).^(1/n-1)).*(-n*abs(y/b).^(n-1)).*sign(y)/b);
    tol = angle-ContactAngle(y);
    if tol > 0        
        y_up = y;
    else
        y_down = y;
    end
    if iter == 100
        break
    end
end

x = (1+p*y).*(1-abs(y/b).^n).^(1/n)*a;

% plot(x,y,'ro','LineWidth',2)
% % quiver(x,y,tp(:,1),tp(:,2),'r')
% x = 0;
% y = 0;