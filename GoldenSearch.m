function ymax = GoldenSearch(a,b,n,p)
p = p/b;
Phi = (-1+sqrt(5))/2;
y0 = [-b b];
f0 = zeros(size(y0));
f = zeros(size(y0));

ymax_old = 0; tol = 1;
while tol>1e-15
    c = y0(2)+Phi*(y0(1)-y0(2));
    y = c;
    f(1) = (1+p*y).*(1-abs(y/b).^n).^(1/n)*a;
    d = y0(1)+Phi*(y0(2)-y0(1));
    y = d;
    f(2) = (1+p*y).*(1-abs(y/b).^n).^(1/n)*a;
%     figure(1);plot(c,f(1),'ro'); hold on; plot(d,f(2),'go');pause
    if f(1)>f(2)
        y0(2) = d;
        f0(2) = f(2);
        c = y0(2)+Phi*(y0(1)-y0(2));
        y = c;
        f(1) = (1+p*y).*(1-abs(y/b).^n).^(1/n)*a;
        d = y0(1)+Phi*(y0(2)-y0(1));
        y = d;
        f(2) = (1+p*y).*(1-abs(y/b).^n).^(1/n)*a;
        ymax_new = c;
    else
        y0(1) = c;
        f0(1) = f(1);
        c = y0(2)+Phi*(y0(1)-y0(2));
        y = c;
        f(1) = (1+p*y).*(1-abs(y/b).^n).^(1/n)*a;
        d = y0(1)+Phi*(y0(2)-y0(1));
        y = d;
        f(2) = (1+p*y).*(1-abs(y/b).^n).^(1/n)*a;
        ymax_new = d;
    end    
    tol = abs(ymax_new-ymax_old);
    ymax_old = ymax_new;
end
ymax = ymax_new;