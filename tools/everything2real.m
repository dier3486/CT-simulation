function A = everything2real(A)
% cast the A to real(A)
%   A = everything2real(A);

classA = class(A);
switch classA
    case 'struct'
        if length(A)>1
            for ii = 1:length(A(:))
                A(ii) = everything2real(A(ii));
            end
        else
            fields = fieldnames(A);
            for ifield = 1:length(fields)
                A.(fields{ifield}) = everything2real(A.(fields{ifield}));
            end
        end
    case 'cell'
        for icell = 1:length(A(:))
            A{icell} = everything2real(A{icell});
        end
    otherwise
        if ~isreal(A)
            A = real(A);
        end
end

end