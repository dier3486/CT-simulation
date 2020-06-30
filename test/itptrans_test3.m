Nx = 512;
xx = (1:Nx) - (Nx+1)/2;
xx_ext = (1:Nx+4) - (Nx+5)/2;

% yy1 = xx+0.3;
Ny = 300;
yy1 = linspace(xx(1), xx(end), Ny) + 0.2;
[index1, alpha1] = interpprepare(xx_ext, yy1);

m = 100;
% X = rand(Nx, 1) + exp(-(xx(:)./Nx.*10).^2);
X = randn(Nx, m);

X_ext = [3*X(1, :)-2*X(2, :); 2*X(1, :)-X(2, :); X; 2*X(end, :)-X(end-1, :); 3*X(end, :)-2*X(end-1, :)];
Y1_line = X_ext(index1(:,1), :).*alpha1(:,1) + X_ext(index1(:,2), :).*alpha1(:,2);

alpha = alpha1(:,2);
beta = 1/2-sqrt(1+4.*alpha.*(1-alpha))./2;
Y1 = X_ext(index1(:,1)-1, :).*(1-alpha+beta)./4 + X_ext(index1(:,1), :).*(2-alpha-beta)./4 + ...
     X_ext(index1(:,2), :).*(1+alpha-beta)./4 + X_ext(index1(:,2), :).*(alpha+beta)./4;

hy = mean(diff(yy1));
k = linspace(0, 0.5, Ny/2+1)./hy;

% b1 = 0.4758;
% b2 = -0.0409;
% 
% A1 = spdiags(repmat([b1/2 1-b1 b1/2], Ny, 1), [-1 0 1], Ny, Ny);
% A2 = spdiags(repmat([b2/2 1-b2 b2/2], Ny, 1), [-1 0 1], Ny, Ny);
% A1(1,1) = 1 - b1/2;
% A1(end,end) = 1 - b1/2;
% A2(1,1) = 1 - b1/2;
% A2(end,end) = 1 - b1/2;
% 
% Z1 = A1\(A2*Y1);
% 
% cut = 0.01;
% k = linspace(0, 0.5, Ny/2+1);
% k = [k, -fliplr(k(2:end-1))];
% f1 = cos(k.*(pi*3)).*(2-sqrt(2))./4 + cos(k.*pi).*(2+sqrt(2))./4;
% f1(f1<cut) = 1;
% f1 = 1./f1(:);
% f1(abs(k)>0.5*hy) = 1;
% 
% Z2 = ifft(fft(Y1).*f1);

