alpha = 0:0.1:1;
k = (0:0.01:0.5)';

h = [1 1 1];
% h = [1 1.5 1];

sigma_a = sqrt((1-alpha).^2 + alpha.^2);
pi2i = 1i*2*pi;
beta1 = 1/2-sqrt(1+4.*alpha.*(1-alpha))./2;
beta2 = 1/8-sqrt(1+8.*alpha.*(1-alpha))./8;
% beta2 = beta2.*0.5;

x1 = exp(-pi2i.*k.*h(1));
x2 = ones(size(k));
x3 = exp(pi2i.*k.*h(2));
x4 = exp(pi2i.*k.*(h(3)+h(2)));

y0 = exp(pi2i.*(k*alpha).*h(2));
y1 = x2*(1-alpha) + x3*alpha;
y2 = (x2 + (x1-x2).*(h(2)/h(1)))*(1-alpha+beta1)./4 + x2*(2-alpha-beta1)./4 + x3*(1+alpha-beta1)./4 + ...
     (x3+(x4-x3).*(h(2)/h(3)))*(alpha+beta1)./4;

y2a = x1*(1-alpha+beta1)./4 + x2*(2-alpha-beta1)./4 + x3*(1+alpha-beta1)./4 + x4*(alpha+beta1)./4;
 
y3 = (x2 + (x1-x2).*(h(2)/h(1)))*beta2 + x2*(1-alpha-beta2) + x3*(alpha-beta2) + (x3+(x4-x3).*(h(2)/h(3)))*beta2;
y3a = x1*beta2 + x2*(1-alpha-beta2) + x3*(alpha-beta2) + x4*beta2;

r1 = y1./y0;
r2 = y2./y0;
r2a = y2a./y0;
r3 = y3./y0;
r3a = y3a./y0;
% h1 = (1-alpha + exp(pi2i.*k)*alpha).*exp(-pi2i.*(k*alpha));


% gamma = -1./sqrt(-alpha.*(1-alpha).*c1+1)./4;
% % gamma = -0.25;
% y4a = y2a + x1*((1-alpha).*gamma) + x2*((alpha.*3-2).*gamma) + x3*((1-alpha.*3).*gamma) + x4*(alpha.*gamma);

c1 = 0.7;
c2 = 1.5;
gamma = c1./sqrt(-alpha.*(1-alpha).*c2+1);
y4 = (x2 + (x1-x2).*(h(2)/h(1)))*((1-alpha).*(1-gamma)+beta1)./4 + x2*(2-alpha-beta1+(2-alpha.*3).*gamma)./4 + ...
      x3*(1+alpha-beta1+(alpha.*3-1).*gamma)./4 + (x3+(x4-x3).*(h(2)/h(3)))*(alpha.*(1-gamma)+beta1)./4;

y4a = x1*((1-alpha).*(1-gamma)+beta1)./4 + x2*(2-alpha-beta1+(2-alpha.*3).*gamma)./4 + ...
      x3*(1+alpha-beta1+(alpha.*3-1).*gamma)./4 + x4*(alpha.*(1-gamma)+beta1)./4;
  


%         t_z_coeff = [ ((1-t_z_alpha(:)).*(1-gamma(:))+beta(:))./4 ...
%                       (2-t_z_alpha(:)-beta(:)+(2-t_z_alpha(:).*3).*gamma(:))./4 ...
%                       (1+t_z_alpha(:)-beta(:)+(t_z_alpha(:).*3-1).*gamma(:))./4 ...
%                       (t_z_alpha(:).*(1-gamma(:))+beta(:))./4];     