alpha = 0:0.05:0.2;
k = (0:0.01:0.5)';

sigma_a = sqrt((1-alpha).^2 + alpha.^2);
pi2i = 1i*2*pi;

h0 = exp(-pi2i.*(k*alpha));
h1 = (1-alpha + exp(pi2i.*k)*alpha).*exp(-pi2i.*(k*alpha));

h2b0 = (1 + sin(2*pi*k)*alpha.*1i).*exp(-pi2i.*(k*alpha));

a_max = 0.2;
h2bm = (1-a_max + a_max.*cos(2*pi*k) + sin(2*pi*k)*alpha.*1i).*exp(-pi2i.*(k*alpha));
b = 0.2;
h2bz = h2bm./(1-b+cos(2*pi*k).*b);
