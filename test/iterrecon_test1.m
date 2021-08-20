imgsize = 256; h = 1;
img0 = phantom('Modified Shepp-Logan',imgsize);
img0 = img0.*5.0;
mu_e = 0.01;
Ne = 1.0e5;
% mu = 0.025;
% Ne = 3.0e5;


Nview = 300;
viewangle = linspace(0, 180, Nview+1);
viewangle = viewangle(1:end-1);
Np = 600;
frequency_scaling = imgsize*sqrt(2)/Np;

P0 = radon(img0, viewangle, Np).*h;
img0b = iradon(P0, viewangle, 'linear', 'hann', frequency_scaling, imgsize);
lenP = size(P0, 1);

I0 = exp(-P0.*mu_e);
I1 = poissrnd(I0.*Ne)./Ne;
P1 = -log(I1)./mu_e;
P1(P1>500) = 500;
img1b = iradon(P1, viewangle, 'linear', 'hann', frequency_scaling, imgsize);

alpha_filt = 1.6;
Kfilt = filterdesign('hann', lenP, 1.0, alpha_filt);
P1f = fconv(P1, Kfilt);
img1a = iradon(P1f, viewangle, 'linear', 'none', frequency_scaling, imgsize);

L = splaplace2D(imgsize);
mu = 5.0;

% F^2
Kf2 = Kfilt.^2.*(pi/2);
% or
Kf0 = Kfilt.*(pi/2);

% iteration
Niter = 50;
alpha = 0.05;
b1 = img1a;
u = zeros(imgsize, imgsize, Niter);
u(:,:,1) = img1b;
rerr = zeros(Niter-1, 1);
for ii = 1:Niter-1
    v1 = radon(u(:,:,ii), viewangle, Np);
    v2 = radon(LLx(L, u(:,:,ii)), viewangle, Np)./mu;
    v = fconv(v1, Kfilt) + fconv(v2, Kf2);
    r = b1 - iradon(v, viewangle, 'linear', 'none', 1.0, imgsize);
    u(:,:,ii+1) = u(:,:,ii) + r.*alpha;
    
    rerr(ii) = sqrt(sum(r(:).^2)).*(1e3/imgsize^2);
end

figure;
plot(log(rerr));


function y = fconv(x, K)

sizeX = size(x);
len = size(K, 1);
x(len, 1) = 0;

y = ifft(fft(x).*K);
y = y(1:sizeX,:);

end

function y = LLx(L, x)

xsize = size(x);
y = reshape(L*x(:), xsize);

end