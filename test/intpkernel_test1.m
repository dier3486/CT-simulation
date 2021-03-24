xx = (-2:0.02:2)';
y1 = zeros(size(xx));
y2 = zeros(size(xx));

c1 = 0.7;
c2 = 1.5;

s = (xx<-1);
alpha = xx(s)+2;
beta = 1/2-sqrt(1+4.*alpha.*(1-alpha))./2;
gamma = c1./sqrt(1-alpha.*(1-alpha).*c2);
y1(s) = (alpha.*(1-gamma)+beta)./4;
y2(s) = (alpha+beta)./4;

s = xx>=-1 & xx<0;
alpha = xx(s)+1;
beta = 1/2-sqrt(1+4.*alpha.*(1-alpha))./2;
gamma = c1./sqrt(1-alpha.*(1-alpha).*c2);
y1(s) = ((1+alpha-beta) + (alpha.*3-1).*gamma)./4;
y2(s) = (1+alpha-beta)./4;

s = xx>=0 & xx<1;
alpha = xx(s);
beta = 1/2-sqrt(1+4.*alpha.*(1-alpha))./2;
gamma = c1./sqrt(1-alpha.*(1-alpha).*c2);
y1(s) = ((2-alpha-beta) + (2-alpha.*3).*gamma)./4;
y2(s) = (2-alpha-beta)./4;

s = xx>=1;
alpha = xx(s)-1;
beta = 1/2-sqrt(1+4.*alpha.*(1-alpha))./2;
gamma = c1./sqrt(1-alpha.*(1-alpha).*c2);
y1(s) = ((1-alpha).*(1-gamma) + beta)./4;
y2(s) = (1-alpha+beta)./4;
