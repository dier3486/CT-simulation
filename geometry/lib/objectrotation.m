function object = objectrotation(object, A)
% rotation by matrix A, A'*A=I.

N = length(object(:));
for ii = 1:N
    object(ii).O = object(ii).O*A;
    object(ii).vector = object(ii).vector*A;
    if isfield(object(ii), 'invV')
        object(ii).invV = A'*object(ii).invV;
    end
end