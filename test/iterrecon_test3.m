imgsize = 512;
img0 = phantom('Modified Shepp-Logan',imgsize);
img0 = img0.*5.0;
img0(img0>2.5) = 2.5;
mu_e = 0.001;
Ndose = 4.0e6;
% mu = 0.025;
% Ne = 3.0e5;

% low density
img0_l = img0(390:430, 210:290);
img0_l(img0_l>1) = 1+(img0_l(img0_l>1)-1).*2e-3.*20.0;
img0(390:430, 210:290) = img0_l;

% img blur
Cblur = 1.5;
img0 = gaussblur(img0, Cblur);

Nview = 500;
viewangle = linspace(0, 180, Nview+1);
viewangle = viewangle(1:end-1);

pbm.delta_d = 0.5;
pbm.h = 0.6;
pbm.Npixel = 700;
pbm.midchannel = (pbm.Npixel+1)/2;
pbm.viewangle = viewangle.*(pi/180);
pbm.Nview = Nview;
pbm.imagesize = imgsize;
pbm.fillmiss = 0;

P0 = parallelprojinimage(pbm, img0, '2D linearinterp');
img0b = filterbackproj2D(P0, pbm, 'hann');

% noise
I0 = exp(-P0.*mu_e);
I1 = poissrnd(I0.*Ndose)./Ndose;
% noise intervene
n = 1024;
c_intv = 0.9;
noise_intv = linspace(1, 0, n/2+1);
noise_intv = [noise_intv fliplr(noise_intv(2:end-1))];
DI = ifft(fft(I1-I0, n, 2).*(1-c_intv*0.8+noise_intv.*c_intv), n, 2);
I2 = I0+DI(:, 1:Nview);

% P
P1 = -log(I1)./mu_e;
P1(P1>500) = 500;
P2 = -log(I2)./mu_e;
P2(P2>500) = 500;

% main filter
alpha_filt = 1.5;
Kfilt = filterdesign('hann', pbm.Npixel, pbm.delta_d, alpha_filt);
img0a = filterbackproj2D(P0, pbm, Kfilt);
img1a = filterbackproj2D(P1, pbm, Kfilt);
img2a = filterbackproj2D(P2, pbm, Kfilt);

% noise filter
alpha_filt = 1.8;
Kfilt_BV = filterdesign('hann', pbm.Npixel, pbm.delta_d, alpha_filt);

% F^2
Kf2 = Kfilt.*Kfilt_BV.*pi;
% or
Kf0 = Kfilt_BV.*pi;
% L
L = splaplace2D(imgsize);
% BV
mu_BV = 0.5;
lambda = 0.03;
mu_F = 3.0;
mu_L = 0.2/mu_F;

% iteration
Niter = 20;
tol_iter = 0.01;
alpha = 0.2;
b1 = img2a;

u = zeros(imgsize, imgsize, Niter);
u(:,:,1) = b1;
u2 = u;
rerr = nan(Niter-1, 1);
fig = figure;
u_BV = [];
for ii = 1:Niter-1
%     u2(:,:,ii) = TVpenalty(u(:,:,ii), mu_BV, lambda);
%     v1 = parallelprojinimage(pbm, u(:,:,ii), '2D linearinterp');
%     v2 = parallelprojinimage(pbm, LLx(L, u(:,:,ii)), '2D linearinterp');
%     v = v1 + fconv(v2, Kf0)./mu;
    % BV
    u_BV = TVpenalty_test1(u(:,:,ii), mu_BV, lambda, u_BV, [800 1200]);
    R_BV = (u(:,:,ii) - u_BV);
    u2(:,:,ii) = u(:,:,ii) - R_BV;
    % style#1 project u; style#2 project u2
    v1 = parallelprojinimage(pbm, u2(:,:,ii), '2D linearinterp');
%     v1 = parallelprojinimage(pbm, u(:,:,ii), '2D linearinterp');
    % L
    R_L = LLx(L, u2(:,:,ii));
    R_BV = R_BV + R_L.*mu_L;
    
    v2 = parallelprojinimage(pbm, R_BV.*mu_F, '2D linearinterp');
    v = v1 + fconv(v2, Kf0);
    r = b1 - filterbackproj2D(v, pbm, Kfilt);
    u(:,:,ii+1) = u(:,:,ii) + r.*alpha;
    
    rerr(ii) = sqrt(sum(r(:).^2)).*(1e3/imgsize^2);
    figure(fig);
    plot(1:Niter-1, rerr,'.-');
    axis([0 Niter 0 max(rerr).*1.1]);
    drawnow;
end
u2(:,:,ii+1) = TVpenalty_test1(u(:,:,ii+1), mu_BV, lambda);

% img2 = u(:,:,ii) + (u2(:,:,ii)-u(:,:,ii)).*tau_BV;





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

