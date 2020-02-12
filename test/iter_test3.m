x0 = phantom.*1000;
theta = 0:0.5:179.5;
emean = 3.0;
[b0, b_noise] = prepare(x0, theta, emean);

x1 = iradon(b0, theta, 'linear', 'Ram-Lak', 256);
x1_noise = iradon(b_noise, theta, 'linear', 'Ram-Lak', 256);

Nv = size(theta(:),1);
Np = size(b0,1);
n = 256^2;
m = Np*Nv;

x_ini = x1_noise;

N = 7;
x = zeros(n, N+1);
x(:, 1) = x_ini(:);
r = zeros(n, N);
a = 0.5;
H1 = zeros(n,1);
H2 = zeros(1,n);
Fkey = 'Ram-Lak';
% Fkey = 'none';
ATb = iradon(b_noise, theta, 'linear', Fkey, 256);
ATb = ATb(:);
alpha = 1.3;
beta =  40.0;
for ii = 1:N
    Ax = radon(reshape(x(:,ii),256,256), theta);
    [FAx, ~] = filterProjections(Ax, Fkey, 1);
    p = panalty1(reshape(x(:,ii),256,256), alpha, beta);
    Ap = radon(p, theta);
    [FAp, ~] = filterProjections(Ap, [Fkey '2'], 1);
    ATAx = iradon(FAx+FAp, theta, 'linear', 'none', 256);
    
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

imtool(reshape(x,256,256*(N+1)),[0, 600]);
imtool(reshape(x-repmat(x1(:),1, N+1),256,256*(N+1)),[-100, 100]);