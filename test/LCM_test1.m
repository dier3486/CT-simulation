% localized-correlation metal artifact fix test1

N = 256;

h = 1.0;

xx = [0:N/2 -N/2+1:-1].*h;
idx = [N/2+2:N 1:N/2+1];

r = 20; ro = 0; rh = 12; 
y0 = 1 - (xx-ro).^2./r^2;
y0(y0<0) = 0;
y0 = sqrt(y0).*rh;

y1 = fft(y0);

w = [0:N/2 -N/2+1:-1]./N./h;

z1 = besselj(1, w.*(r*pi*2))./w.*(rh/2).*exp(-(1.0i*2*pi*ro).*w);
z1(1) = r*rh*pi/2;

z2 = besselj(1, w.*(r*pi*2)).*sign(w).*(rh/2).*exp(-(1.0i*2*pi*ro).*w);