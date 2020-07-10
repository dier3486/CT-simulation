% img0=load("F:\data-Dier.Z\imgdata\Zfilter\ImageOut.mat");

N = 64;
x = -N/2+1:N/2;
t = [0:N/2 -N/2+1:-1];

y0 = zeros(1, N);
% y0(30:38) = 50;

% y0 = y0 + randn(1, N);
y0(1:32) = img0.ImageOut(52, 280, :)-1000;

a1 = 2;
a2 = 1.0;
z1=exp(-x.^2./a1).*sqrt(1/a1/pi);
z2=exp(-x.^2./a2).*sqrt(1/a2/pi);

% wf1 = 28;
% f1 = exp(-t.^2./wf1^2);

w2 = 0.025;
f1 = sinc(t.*w2);

% w2 = 1.0;
% f1 = ifft((-1).^t.*exp(-abs(t).*w2));
% f1 = real(f1(1)./f1);

% f1 = real(fft(z0)./fft(y0));
% f1 = real(f1(1)./f1);

% f1 = real(fft(z2)./fft(z1));
% f1 = f1(1)./f1;

f2 = ifft(fft(f1).^2);
% I1 = ones(1, N);
lambda = 0.0;

k = (f1+lambda)./(f1.^2+lambda);
y1 = real(ifft(fft(y0).*k));

A1 = reshape(img0.ImageOut, [], 32);
A1 = [A1 fliplr(A1)];
A2 = real(ifft(fft(A1, [], 2).*k, [], 2));