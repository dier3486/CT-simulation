function [index, alpha] = interpprepare(X, Xq, extrap)
% prepare the interpolation
% [index, alpha] = interpprepare(X, Xq, extrap);

if nargin<3
    extrap = 0;
end

N = length(X);
Nq = length(Xq);
index = zeros(2, Nq);
alpha = zeros(2, Nq);

index0 = sum(X(:)' < Xq(:), 2);
index0(Xq==X(1)) = 1;
index(1,:) = index0;
index(2,:) = index0+1;

s = (index0>0) & (index0<N);
alpha(2, s) = (Xq(s) - X(index0(s)))./(X(index0(s)+1)- X(index0(s)));
alpha(1, s) = 1-alpha(2, s);

switch extrap
    case 0
        index(:, ~s) = 1;
    case 'nan'
        index(:, ~s) = 1;
        alpha(:, ~s) = nan;
    case 'extrap'
        index(1, index0==0) = 1;
        index(2, index0==0) = 2;
        index(1, index0==N) = N-1;
        index(2, index0==N) = N;
        alpha(2, ~s) = (Xq(~s) - X(index(1, ~s)))./(X(index(2, ~s))- X(index(1, ~s)));        
        alpha(1, ~s) = 1-alpha(2, ~s);
    otherwise
        index(:, ~s) = 1;
        alpha(:, ~s) = nan;
end
