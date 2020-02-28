function CG = crsmatrix_test1(p)

p = p(:);
N = size(p,1);

P = diag(p);
D = diag(ones(N,1)) -diag(ones(N-1,1),1);

A = eye(N) - D*P*D';
Ad = diag(diag(A));
Ap = A-Ad;

C = (eye(N)-Ap)*inv(Ad);

g = C\ones(N,1);
CG = C*diag(g);

% % before gain norm
% a = 0.1;
% m = 6;
% n = 12;
% p = ones(n-1, 1).*a;
% p(m:m:end) = 0;
% 
% Ad = diag(ones(n,1)-a.*2);
% Ap = diag(p, -1)+diag(p, 1);
% A = Ad+Ap;
% 
% C = (eye(N)-Ap)*inv(Ad) - eye(N);
