function A = everything2single(A, fromwhat, towhat)
% cast the A's values in class fromwhat to class towhat
% A = everything2single(A, fromwhat, towhat)
% or A = everything2single(A) to cast all double floats to single floats

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
        A = cast(A, towhat);
    otherwise
        if strcmp(fromwhat, 'everything')
            A = cast(A, towhat);
        end
end

return