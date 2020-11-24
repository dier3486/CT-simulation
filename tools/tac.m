function varargout = tac(A, dims)
% inverse cat
% varargout = tac(A, dims);
% or varargout = tac(A);
% e.g. [a b c] = tac([1 2 3]); to return a=1 b=2 c=3.
% the 'dims' is just like what input in num2cell.
% [a1 a2] = tac(A, [1 2]); to return a1=A(:,:,1), a2=A(:,:,2);

if nargin>1
    Ac = num2cell(A, dims);
else
    Ac = num2cell(A);
end
Ac = Ac(:);

if nargout > length(Ac)
    Ac{nargout} = [];
end

varargout = Ac;

end
