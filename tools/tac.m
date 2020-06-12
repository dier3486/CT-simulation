function varargout = tac(A, dims)
% inverse cat
% varargout = tac(A, dims);
% or varargout = tac(A);
% e.g. [a b c] = tac([1 2 3]); to return a=1 b=2 c=3.

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
