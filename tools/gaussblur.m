function B = gaussblur(A, sigma)
% Gauss blur
% B = gaussblur(A, sigma);

% complex recurse
if ~isreal(sigma)
    B = complex(gaussblur(real(A), real(sigma)),  gaussblur(imag(A), imag(sigma)));
    return
end
% Well, only 'real' sigma will be used in blur.

Asize = size(A);

if length(sigma(:)) == 1
    sigma = [sigma sigma];
end
% no 0 simga
if all(abs(sigma)<eps)
    B = A;
    return;
end

m = 2.^nextpow2(Asize(1));
n = 2.^nextpow2(Asize(2));

K = gausskernal2(m, n, sigma);

if isreal(A)
    symflag = 'symmetric';
else
    symflag = 'nonsymmetric';
end

B = ifft2(fft2(A, m, n).*K, symflag);
B = B(1:Asize(1), 1:Asize(2), :);

end


function K = gausskernal2(m, n, sigma)

xx = [0:m/2 fliplr(-1:-1:-m/2+1)];
yy = [0:n/2 fliplr(-1:-1:-n/2+1)];
K = real(fft2(exp(-xx(:).^2./(2*sigma(1))-yy.^2./(2*sigma(2))))./(pi*2*sqrt(sigma(1)*sigma(2))));
K = K./K(1,1);

end