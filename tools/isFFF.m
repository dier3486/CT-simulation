function r = isFFF(X)
% check if X is all F, 0xFF..FF

r = reshape(all(reshape(typecast(X(:), 'uint8')==255, [], length(X(:)))), size(X));