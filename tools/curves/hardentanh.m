function y = hardentanh(x, sigma, c)
% harden tanh
% y = hardentanh(x, sigma, c).
%   0<sigma<=1
%   and hardentanh(x, 1, 1) == tanh(x);

if nargin<3
    c = 1.0;
end

if abs(1-sigma) < eps
    % it is tanh
    y = (1-exp(-x.*2./c))./(1+exp(-x.*2./c)).*c;
else
    b = (1+sigma)/(1-sigma)*log(sigma)/c;
    c = -c./log(sigma); 
    y = log((1 + sigma.*exp(x.*b))./(sigma + exp(x.*b))).*c;
end

end