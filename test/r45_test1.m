% m to n

m = 363;
n = 512+1;

xygrid = (1:m) - (m+1)/2;
[Xm, Ym] = ndgrid(xygrid);
XY = [Xm(:) Ym(:)]*[1 1; -1 1] + (n+1)/2;


sXY = any((XY<1) | (XY>n), 2);
S = (XY(:,2)-1).*n + XY(:,1);
S(sXY) = n^2+1;

% A = magic(m);
% B1 = zeros(n^2+1,1);
% B1(S) = A;
% B = reshape(B1(1:end-1), n, n);


img_0 = reshape(img_363{1}, m^2, []);
Nimg = size(img_0, 2);
img_1 = zeros(n^2+1, Nimg, 'single');
img_1(S, :) = img_0;
img_1 = reshape(img_1(1:end-1, :), n, n, Nimg);

C = [0 1 0; 1 0 1; 0 1 0]./4;
for ii = 1:Nimg
    img_tmp = img_1(:,:,ii);
    img2 = conv2(img_tmp, C, 'same');
    img_tmp(2:2:end) = img2(2:2:end);
    img_1(:,:,ii) = img_tmp;
end
