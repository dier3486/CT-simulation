function y = TVexp(x, a)

% y = (1-exp(-x./a))./(1+exp(-x./a));
% y(x./a>10) = -1;
% y(x./a<-10) = 1;

% y = (abs(x).^(2*a) + abs(x).^a.*(1+a))./(abs(x).^a+1).^2;
% y = y.*sign(x);

y = abs(x).^(a-1)./a.*sign(x);
