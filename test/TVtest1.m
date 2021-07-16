% img0 = loaddicomimg('F:\Data\extra\Insitum528BrainImage\Insitum528_brain_tilt_1.25mm\ImageData_bioON_brainPlus\2.1620612318094.pd.2.1620612318096');
% img0b = loaddicomimg('F:\Data\extra\Insitum528BrainImage\Insitum528_brain_tilt_1.25mm\ImageData_bioON_soft1\2.1620612318094.pd.2.1620612318096');

f0 = img0(:,:,55);
f0b = img0b(:,:,55);
imgsize = size(f0);

Crange = 1045+[-100 100];
s1 = (f0>=Crange(1)) & (f0<=Crange(2));
f1 = f0;
f1(f1<Crange(1)) = Crange(1);
f1(f1>Crange(2)) = Crange(2);

Cblur = 5.0;
Cgw = 1035;

f2 = gaussblur(f1, Cblur);
s2 = f2>1035;

f3 = f1.*s2 + Crange(1).*(~s2);

Niter = 30;
lambda = 0.1;
mu1 = 0.7;
mu2 = 0.4;
mu = mu2.*s2 + mu1.*(~s2);

% u = zeros([imgsize Niter]);
b = zeros([imgsize 2]);
d = zeros([imgsize 2]);
delta = zeros(1, Niter);

for ii = 1:Niter
    if ii == 1 
        u0 = f1;
    else
        u0 = u;
    end
    u = funG(f1, u0, b, d, mu, lambda);
    delta(ii) = sqrt(sum((u(:)-u0(:)).^2));
    [d, b] = fundbyub(u, b, lambda);
end


function [d1, b1] = fundbyub(u, b0, lambda)
% d_x^{k+1} = shrink(D_xu^{k+1}+b_x^{k}, 1/\lambda)
% d_y^{k+1} = shrink(D_yu^{k+1}+b_y^{k}, 1/\lambda)

[nx, ny] = size(u);

du = zeros(nx, ny, 2);
du(:,:,1) = [u(2:end,:) - u(1:end-1, :); zeros(1, ny)];
du(:,:,2) = [u(:, 2:end) - u(:, 1:end-1), zeros(nx, 1)];

b0 = b0 + du;
d1 = max(abs(b0)-1/lambda, 0).*(b0./abs(b0));
d1 = fillmissing(d1, 'constant', 0);

b1 = b0 - d1;

end

function u1 = funG(f, u, b, d, mu, lambda)
[nx, ny] = size(u);

u1 = u([2:nx nx], :) + u([1 1:nx-1], :) + u(:, [2:ny ny]) + u(:, [1 1:ny-1]);
u1 = u1 + d([1 1:nx-1], :, 1) - d(:, :, 1) + d(:, [1 1:ny-1], 2) - d(:,:,2);
u1 = u1 - b([1 1:nx-1], :, 1) + b(:, :, 1) - b(:, [1 1:ny-1], 2) + b(:,:,2);
u1 = u1.*(lambda./(mu+4*lambda)) + f.*(mu./(mu+4*lambda));

end