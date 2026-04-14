function x = complex2float(x)
% trans x to (real(x); imag(x))

xsize = size(x);
xsize(1) = xsize(1)*2;

x = reshape([real(x(:)) imag(x(:))]', xsize);