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
K = gausskernal(Asize, sigma);

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

xx = [0:Asize(1)/2 fliplr(-1:-1:-Asize(1)/2+1)].*(sigma(1)/Asize(1));
yy = [0:Asize(2)/2 fliplr(-1:-1:-Asize(2)/2+1)].*(sigma(2)/Asize(2));
K = exp(-(xx(:).^2+yy.^2).*(pi^2/2));

end