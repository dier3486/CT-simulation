function B = gaussblur(A, sigma)
% Gauss blur
% B = gaussblur(A, sigma);

Asize = size(A);
Ndim = ndims(A);
if Ndim==2 && any(Asize==1)
    Ndim = 1;
end
if length(sigma(:)) == 1
    sigma = [sigma sigma];
end
% no 0 simga
if all(abs(sigma)<eps)
    B = A;
    return;
end

% K = gausskernal(Asize, sigma);
K = gausskernal2(Asize, sigma);

switch Ndim
    case 1
        B = ifft(fft(A).*K, 'symmetric');
    case 2
        B = ifft2(fft2(A).*K, 'symmetric');
    otherwise
        A = reshape(A, Asize(1), Asize(2), []);
        B = zeros(size(A));
        for ii = 1:size(A,3)
            B(:,:,ii) = ifft2(fft2(A(:,:,ii)).*K, 'symmetric');
        end
end
B = reshape(B, Asize);        

end


function K = gausskernal(Asize, sigma)

xx = [0:Asize(1)/2 fliplr(-1:-1:-Asize(1)/2+1)]./Asize(1);
yy = [0:Asize(2)/2 fliplr(-1:-1:-Asize(2)/2+1)]./Asize(2);
K = exp(-(xx(:).^2.*sigma(1) + yy.^2.*sigma(2)).*(pi^2*2));

end


function K = gausskernal2(Asize, sigma)

xx = [0:Asize(1)/2 fliplr(-1:-1:-Asize(1)/2+1)];
yy = [0:Asize(2)/2 fliplr(-1:-1:-Asize(2)/2+1)];
K = real(fft2(exp(-xx(:).^2./(2*sigma(1))-yy.^2./(2*sigma(2))))./(pi*2*sqrt(sigma(1)*sigma(2))));
K = K./K(1,1);

end