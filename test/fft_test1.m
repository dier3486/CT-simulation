dx = 0.1;
Fs = 1/dx;
Nx = 40;
xx = (-Nx:Nx).*dx;

sigma = 0.5;
y1 = (1/sqrt(2*pi)/sigma).*exp(-xx.^2./(sigma^2*2));

N = 2^ceil(log2(Nx*2+1))*2;
y2 = zeros(1, N);
y2(1:Nx+1) = y1(Nx+1:end);  y2(N-Nx+1:end) = y1(1:Nx);

fy2 = fft(y2).*dx;
f = Fs*(0:(N/2))/N;

fy0 = exp(-sigma^2.*(f.*pi.*2).^2./2);



