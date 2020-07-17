alpha = 0:0.1:1;
k = (0:0.01:0.5)';

sigma_a = sqrt((1-alpha).^2 + alpha.^2);
pi2i = 1i*2*pi;

h0 = exp(-pi2i.*(k*alpha));
h1 = (1-alpha + exp(pi2i.*k)*alpha).*exp(-pi2i.*(k*alpha));

beta = 1/2-sqrt(1+4.*alpha.*(1-alpha))./2;
h2 = (exp(-pi2i.*k)*(1-alpha+beta).*(1/4) + (1-alpha+1-beta)./4 + exp(pi2i.*k)*(alpha+1-beta)./4 + ...
     exp(pi2i.*2.*k)*(alpha+beta).*(1/4)).*exp(-pi2i.*(k*alpha));


gamma3 = 0.55;
%  h3 = (exp(-pi2i.*2.*k)*(-1+alpha-beta).*(gamma/8) + ...
%       exp(-pi2i.*k)*(2-alpha.*(2+gamma)+beta.*(2+3*gamma))./8 + ...
%       (2+gamma - (alpha+beta).*(1+gamma))./4 + ...
%       exp(pi2i.*k)*(1 + (alpha-beta).*(1+gamma))./4 + ...
%       exp(pi2i.*2.*k)*(alpha.*(2+gamma) + beta.*(2+gamma*3) - gamma)./8 + ...
%       exp(pi2i.*3.*k)*(-alpha-beta).*(gamma/8)).*exp(-pi2i.*(k*alpha));

h3 = exp(-pi2i.*2.*k) * (-(1-alpha+beta).*(gamma3/8)) + ...                              % -2
     exp(-pi2i.*k) * ((1-alpha+beta)./4 + (-alpha+beta.*3).*(gamma3/8)) + ...         % -1
     1 .* ((2-alpha-beta)./4 + (1-alpha-beta).*(gamma3/4)) + ...             % 0
     exp(pi2i.*k) * ((1+alpha-beta)./4 + (alpha-beta).*(gamma3/4)) + ...               % +1
     exp(pi2i.*2.*k) * ((alpha+beta)./4 + (alpha-1+beta.*3).*(gamma3/8)) + ...           % +2
     exp(pi2i.*3.*k) * (-(alpha+beta).*(gamma3/8));
h3 = h3.*exp(-pi2i.*(k*alpha));

% g2 = fzero(@(x) (g1/8)^2+(1/4-x/4)^2+(1/2+g1/4+x*3/8)^2+1/16+(g1/8+x/8)^2-1, 0.5);
gamma4 = 0.65;
h4 = (exp(-pi2i.*k) * ((1-alpha+beta)./4 + (-2+alpha+beta).*(gamma4/8)) + ...
      1 .* ((1-alpha+1-beta)./4 + (3-alpha.*3-beta).*(gamma4/8))+ ...
      exp(pi2i.*k) * ((alpha+1-beta)./4 + (alpha.*3-beta).*(gamma4/8))+ ...
      exp(pi2i.*2.*k) * ((alpha+beta)./4 + (-1-alpha+beta).*(gamma4/8))) ...
      .*exp(-pi2i.*(k*alpha));
  
g1 = (sqrt(19)-2)*2/3;
g2 = (sqrt(170)-4)/7;
gc = 0.5;
gs = 0.43;
gamma5a = g1*gc*gs;
gamma5b = g2*(1-gc)*gs;

h5 = exp(-pi2i.*2.*k) * (-(1-alpha+beta).*(gamma5a/8)) + ...                              % -2
     exp(-pi2i.*k) * ((1-alpha+beta)./4 + (-alpha+beta.*3).*(gamma5a/8) + (-2+alpha+beta).*(gamma5b/8)) + ...         % -1
     1 .* ((2-alpha-beta)./4 + (1-alpha-beta).*(gamma5a/4) + (3-alpha.*3-beta).*(gamma5b/8)) + ...             % 0
     exp(pi2i.*k) * ((1+alpha-beta)./4 + (alpha-beta).*(gamma5a/4) + (alpha.*3-beta).*(gamma5b/8)) + ...               % +1
     exp(pi2i.*2.*k) * ((alpha+beta)./4 + (alpha-1+beta.*3).*(gamma5a/8) + (-1-alpha+beta).*(gamma5b/8)) + ...           % +2
     exp(pi2i.*3.*k) * (-(alpha+beta).*(gamma5a/8));
h5 = h5.*exp(-pi2i.*(k*alpha));