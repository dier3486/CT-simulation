A = squeeze(df0.rawdata(:, 27, :));

% h = 0.02;
% % L = 0.8;
% a = h.*0.5;
% xx = [0:h:0.5 (-0.5+h):h:-h];
% N = length(xx);
% 
% 
% y1 = exp(-xx.^2./a).*sqrt(1/a/pi);
% y2 = -xx.*exp(-xx.^2./a).*sqrt(1/a^3/pi).*2;
% 
% fy1 = fft(y1).*(1/N);
% fy2 = fft(y2).*(1/N);
% 
% zz = xx.*N;
% u1 = exp(-(pi.*zz).^2.*a);
% u2 = zz.*exp(-(pi.*zz).^2.*a).*(pi*2);

Np = size(A, 1);
N = 1024;
% Nv = 1152*2;
Nv = 2304;
f0 = A(:,1:Nv).*(pi/2);
f1 = fft(f0, N);

h = 1/N;
xx = linspace(0, 1, N+1);
xx = [xx(1:N/2+1) xx(N/2+2:N)-1];

a = (h*1.0)^2;
y1 = exp(-xx.^2./a).*sqrt(1/a/pi);
y2 = -xx.*exp(-xx.^2./a).*sqrt(1/a^3/pi).*2;

fy1 = fft(y1).*(1/N);
fy2 = fft(y2).*(1/N);

zz = xx.*N;
u1 = exp(-(pi.*zz).^2.*a);
u2 = zz.*exp(-(pi.*zz).^2.*a).*(pi*2);

uu2 = ifft(-u2.^2);
f2 = ifft(f1.*(1i.*u2(:)));

f3 = real([zeros(1, Nv); f2(2:end-1, :) - f2(1:end-2, :)./2 - f2(3:end, :)./2; zeros(1, Nv)])./uu2(1);
% f3m = median(f3, 2);
f3m = median(f3, 2)./std(f3, [], 2)./3;
% f3m = mean(f3, 2)./std(f3, [], 2)./3;
f3m(Np+1:end) = 0;
f3m(1) = 0;

f4 = ifft(fft(f3m).*(1i.*u2(:)));

g1 = (f0 - real(f4(1:Np)))./(pi/2);
