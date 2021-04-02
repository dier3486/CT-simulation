function [y, f] = whatisfft(dx, lubond, fun, varargin)
% A sample of how to use fft. Anyone puzzled by fft could get help from this ;)
% [y, f] = whatisfft(dx, lubond, fun, vargin);
% e.g. [y, f] = whatisfft(0.1, [-10 10], @(x, a, b) exp(-(x-b).^2./(2*a^2))./(sqrt(2*pi)*a) , 1.1, 2); to get the fft
% transformation of a normal distribution curve. 
% And plot(f, abs(y)) to show the curve, where the f is in Hz, and note omiga=f*(2*pi) 

% xx
xx = lubond(1):dx:lubond(2);
Nx = length(xx);
N = 2^ceil(log2(Nx));
% fx = fun(xx)
if isa(fun, 'function_handle')
    fx = fun(xx, varargin{:});
else
    fx = fun;
end
% fft
y = fft(fx, N);
% norm
y = y(1:N/2+1).*dx;
f = (0:(N/2))/N./dx;
end