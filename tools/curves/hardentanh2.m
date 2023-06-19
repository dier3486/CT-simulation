function y = hardentanh2(x, sigma, c)
% harden tanh, stable and vector (c) patch of hardentanh().
% y = hardentanh2(x, sigma, c).
%   0<sigma<=1
%   and hardentanh(x, 1, 1) == tanh(x);

if nargin<3
    c = x.*0 + 1;
end

if abs(1-sigma) < eps
    % it is tanh
    y = tanh(x./c).*c;
else
    logsigma = log(sigma);
    xb = x./c .* ((1+sigma)/(1-sigma)*logsigma);
    a = 1/tanh(logsigma);
    y = log(sigma/(a + 1).*(tanh((xb-logsigma)./2) + a)).*(-c./logsigma);
end

% I know y is nan only when x==0 & c==0.
y = fillmissing(y, 'constant', 0);

end