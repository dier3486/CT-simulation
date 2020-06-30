alpha = 0:0.1:1;
k = (0:0.01:1.0)';

sigma_a = sqrt((1-alpha).^2 + alpha.^2);
pi2i = 1i*2*pi;

h0 = exp(-pi2i.*(k*alpha));
h1 = (1-alpha + exp(pi2i.*k)*alpha).*exp(-pi2i.*(k*alpha));

beta = 1/2-sqrt(1+4.*alpha.*(1-alpha))./2;
h2 = (exp(-pi2i.*k)*(1-alpha+beta).*(1/4) + (1-alpha+1-beta)./4 + exp(pi2i.*k)*(alpha+1-beta)./4 + ...
     exp(pi2i.*2.*k)*(alpha+beta).*(1/4)).*exp(-pi2i.*(k*alpha));

% % b = max(beta);
% b = 1/3;
% h2z = h2./(1-b+cos(2*pi*k).*b);

% a = (1-alpha).^2 + alpha.^2;
% b1 = (a-1)/2;
% b2 = (a-1)/2;
% 
% h3 = (exp(-pi2i.*k)*(1-alpha+b1).*(1/4) + (2-alpha-2.*b1+b2)./4 + exp(pi2i.*k)*(1+alpha+b1-b2.*2)./4 + ...
%      exp(pi2i.*2.*k)*(alpha+b2).*(1/4)).*exp(-pi2i.*(k*alpha));
 
% bz = max(-beta)*2
% bz1 = 0.45;
bz1 = (2-sqrt(2))/3;
h2z1 = h2./(1-bz1+cos(2*pi*k).*bz1);

f1 = cos(k.*(pi*3)).*(2-sqrt(2))./4 + cos(k.*pi).*(2+sqrt(2))./4;
% f1 = abs(h2(:, 1));

bz2 = lsqnonlin(@(x) 1-sum(x)+x(1).*cos(k.*(pi*2))+x(2).*cos(k.*(pi*4)) - f1,[0 0]);
f2 = 1-sum(bz2)+bz2(1).*cos(k.*(pi*2))+bz2(2).*cos(k.*(pi*4));

indexk = 1:51;

bz3 = lsqnonlin(@(x) (1-x(1)+x(1).*cos(k(indexk).*(pi*2)))./(1-x(2)+x(2).*cos(k(indexk).*(pi*2))) - f1(indexk),[0 0]);
f3 = (1-bz3(1)+bz3(1).*cos(k.*(pi*2)))./(1-bz3(2)+bz3(2).*cos(k.*(pi*2)));

% ff1 = real(fft([f1; flipud(f1(2:end-1))]));

% hf1 = bz3(1);
hf1 = 0.5;

A2=full(spdiags(repmat([hf1/2 1-hf1 hf1/2], 16, 1), [-1 0 1], 16, 16));
A2(1,1)=1-hf1/2;
A2(end,end)=1-hf1/2;