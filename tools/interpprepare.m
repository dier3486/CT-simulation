function [index, alpha] = interpprepare(X, Xq, extrap)
% prepare the interpolation
% [index, alpha] = interpprepare(X, Xq, extrap);
% then the Vq = interp1(X,V,Xq,'linear',extrap) is equavelant to 
% Vq = V(index(1,:)).*alpha(1,:) + V(index(2,:)).*alpha(2,:);

if nargin<3
    extrap = 0;
end

N = length(X);
Nq = length(Xq);
index = zeros(Nq, 2);
alpha = zeros(Nq, 2);

index0 = sum(X(:)' < Xq(:), 2);
index0(Xq==X(1)) = 1;
index(:, 1) = index0;
index(:, 2) = index0+1;

s = (index0>0) & (index0<N);
alpha(s, 2) = (Xq(s) - X(index0(s)))./(X(index0(s)+1)- X(index0(s)));
alpha(s, 1) = 1-alpha(s, 2);

switch extrap
    case 0
        index(:, ~s) = 1;
    case 'nan'
        index(:, ~s) = 1;
        alpha(:, ~s) = nan;
    case 'extrap'
        index(index0==0, 1) = 1;
        index(index0==0, 2) = 2;
        index(index0==N, 1) = N-1;
        index(index0==N, 2) = N;
        alpha(~s, 2) = (Xq(~s) - X(index(~s, 1)))./(X(index(~s, 2))- X(index(~s, 1)));        
        alpha(~s, 1) = 1-alpha(~s, 2);
    otherwise
        index(~s, :) = 1;
        alpha(~s, :) = nan;
end
