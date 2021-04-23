
Nx = 1000;
% x0 = rand(Nx,1);
x0 = linspace(0.1,1,Nx)';
a0 = [0.7 -0.8 1.0];
Na = length(a0);
y0 = polyval([a0 0], x0);
x1 = x0 + randn(Nx,1).*0.1;
p1 = polyfit(x1, y0, Na);
% p1 = [polyfit(x1, y0./x1, Na-1) 0];

A = x1.^(Na:-1:1);
p2 = A\y0;
p2 = [p2' 0];

u_ini = zeros(1, Na);
u_ini(end) = 1;
Niter = 10;
% u = zeros(Niter+1, Na);
% u(1, :) = u_ini;
u = u_ini;
x_ii = [];
alpha = 1.0;
tol = 1e-8;
for iter = 1:Niter
%     u_ii = u(iter, :);
    x_ii = invpolyval([u 0], y0, x_ii);
    r = x_ii - x1;
    
    du = (Na:-1:1).*u;
    dy = polyval(du, x_ii);
    A = (x_ii.^(Na:-1:1))./dy;
    AA = A'*A;
    r_u = AA\(A'*r);
    u = u + r_u'.*alpha;
    if all(abs(r_u)<tol)
        break;
    end
end

xx = linspace(0,1,100);
y1 = polyval([a0 0], xx);
y2 = polyval(p2, xx);
y3 = polyval([u 0], xx);

figure; hold on
plot(x1, y0, 'c.', 'MarkerSize', 2.0);
plot(xx, y1, 'g--', 'LineWidth', 2.0);
plot(xx, y2, 'r', 'LineWidth', 1.0);
plot(xx, y3, 'b',  'LineWidth', 1.0);
axis equal
grid on