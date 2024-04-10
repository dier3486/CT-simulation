function y = hardentanh(x, sigma, c, sigmalog_flag)
% harden tanh, robust version
% y = hardentanh(x, sigma, c).
%   0<sigma<=1
%   and hardentanh(x, 1, 1) == tanh(x);

if nargin<3
    c = ones(1,2, 'like', x);
end
if length(c) == 1
    c = [c c];
end
% to set the sigma by log(sigma)
if nargin<4
    sigmalog_flag = false;
end

if ~sigmalog_flag
    logsigma = log(sigma);
else
    logsigma = sigma;
    sigma = exp(logsigma);
end

% c = |c|
c = abs(c);

spos = x > 0;
sneg = x < 0;
y = x;

if abs(1-sigma) < eps(classGPU(sigma))*100
    % it is almost tanh
    b = -(1+sigma)./c;
    y(spos) = (1 - exp(x(spos).*b(2))) ./ (sigma + exp(x(spos).*b(2))) .* c(2);
    y(sneg) = (exp(-x(sneg).*b(1)) - 1) ./ (sigma.*exp(-x(sneg).*b(1)) + 1) .* c(1);

elseif sigma < eps(classGPU(sigma))
    % deap small sigma, almost zig-line
    b = (1+sigma)/(1-sigma)*logsigma./c;
    y(spos) = (log(1 + sigma.*exp(x(spos).*b(2))) - log(sigma + exp(x(spos).*b(2)))) .* (-c(2)/logsigma);
    y(sneg) = (log(exp(-x(sneg).*b(1)) + sigma) - log(sigma.*exp(-x(sneg).*b(1)) + 1)) .* (-c(1)/logsigma);

else
    % normal case
    b = (1+sigma)/(1-sigma)*logsigma./c;
    y(spos) = log( (1 + sigma.*exp(x(spos).*b(2))) ./ (sigma + exp(x(spos).*b(2))) ) .* (-c(2)/logsigma);
    y(sneg) = log( (exp(-x(sneg).*b(1)) + sigma) ./ (sigma.*exp(-x(sneg).*b(1)) + 1) ) .* (-c(1)/logsigma);
end

end