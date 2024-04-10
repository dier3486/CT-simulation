function y = hardentanh2(x, sigma, c, sigmalog_flag)
% harden tanh, vector c patch of hardentanh().
% y = hardentanh2(x, sigma, c).
%   0<sigma<=1
%   and hardentanh(x, 1, 1) == tanh(x);

if nargin<3
    c = x.*0 + 1;
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
    y(spos) = (1 - exp(x(spos).*b(spos))) ./ (sigma + exp(x(spos).*b(spos))) .* c(spos);
    y(sneg) = (exp(-x(sneg).*b(sneg)) - 1) ./ (sigma.*exp(-x(sneg).*b(sneg)) + 1) .* c(sneg);

elseif sigma < eps(classGPU(sigma))
    % deap small sigma
    b = (1+sigma)/(1-sigma)*logsigma./c;
    y(spos) = (log(1 + sigma.*exp(x(spos).*b(spos))) - log(sigma + exp(x(spos).*b(spos)))) .* (-c(spos)/logsigma);
    y(sneg) = (log(exp(-x(sneg).*b(sneg)) + sigma) - log(sigma.*exp(-x(sneg).*b(sneg)) + 1)) .* (-c(sneg)/logsigma);

else
    % normal case
    b = (1+sigma)/(1-sigma)*logsigma./c;
    y(spos) = log( (1 + sigma.*exp(x(spos).*b(spos))) ./ (sigma + exp(x(spos).*b(spos))) ) .* (-c(spos)/logsigma);
    y(sneg) = log( (exp(-x(sneg).*b(sneg)) + sigma) ./ (sigma.*exp(-x(sneg).*b(sneg)) + 1) ) .* (-c(sneg)/logsigma);
end

end