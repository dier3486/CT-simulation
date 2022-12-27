function A = everything2single(A, fromwhat, towhat)
% cast the A's values in class fromwhat to class towhat
% A = everything2single(A, fromwhat, towhat)
% or A = everything2single(A) to cast all double floats to single floats
% which equivalent to  A = everything2single(A, 'double', 'single').
% To use A = everything2single(A, 'any', 'double') to cast all numerical
% value to double

if nargin<2
    fromwhat = 'double';
elseif isempty(fromwhat)
    fromwhat = 'double';
end

if nargin<3
    towhat = 'single';
elseif isempty(towhat)
    towhat = 'single';
end

classA = class(A);
switch classA
    case 'struct'
        if length(A)>1
            for ii = 1:length(A)
                A(ii) = everything2single(A(ii), fromwhat, towhat);
            end
        else
            fields = fieldnames(A);
            for ifield = 1:length(fields)
                A.(fields{ifield}) = everything2single(A.(fields{ifield}), fromwhat, towhat);
            end
        end
    case 'cell'
        for icell = 1:length(A)
            A{icell} = everything2single(A{icell}, fromwhat, towhat);
        end
    case fromwhat
        A = cast2what(A, towhat);
    otherwise
        if strcmpi(fromwhat, 'any') && isnumeric(A)
            A = cast2what(A, towhat);
        end
end

end

function A = cast2what(A, towhat)

switch lower(towhat)
    case 'gpuarray'
        A = gpuArray(A);
    case 'gpusingle'
        A = gpuArray(single(A));
    case 'gather'
        A = gather(A);
    otherwise
        A = cast(A, towhat);
end
end