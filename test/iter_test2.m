x0 = phantom.*1000;
theta = 0:0.5:179.5;
b = radon(x0, theta);

Nv = size(theta(:),1);
Np = size(b,1);
n = 256^2;
m = Np*Nv;

x_ini = iradon(b, theta, 'linear', 'Ram-Lak', 256);

N = 10;
x = zeros(n, N+1);
x(:, 1) = x_ini(:);
r = zeros(n, N);
a = 1.0;
H1 = zeros(n,1);
H2 = zeros(1,n);
Fkey = 'Ram-Lak';
% Fkey = 'none';
ATb = iradon(b, theta, 'linear', Fkey, 256);
ATb = ATb(:);
for ii = 1:N
    Ax = radon(reshape(x(:,ii),256,256), theta);
    ATAx = iradon(Ax, theta, 'linear', Fkey, 256);
    r(:,ii) =  ATAx(:) - ATb;
    if ii>1
        y = r(:,ii) - r(:,ii-1);
        Hy = y + H1*(H2*y);
        V1 = (dx-Hy)./(dx'*Hy);
        H1 = [H1 + V1*(dx'*H1), V1];
        H2 = [H2; dx'];
        dx = -(r(:, ii) + H1*(H2*r(:,ii))).*a;
    else
        dx = -r(:, ii).*a;
    end
    
    x(:,ii+1) = x(:,ii) + dx;
end

err = sqrt(sum(r.^2,1));
figure;
plot(err, '.-');