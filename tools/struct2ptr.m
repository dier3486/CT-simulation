function A = struct2ptr(A)

fields = fieldnames(A);
for ii = 1:length(fields)
    fclass = class(A.(fields{ii}));
    switch fclass
        case 'struct'
           A.(fields{ii}) =  struct2ptr(A.(fields{ii}));
        case 'char'
            if length(A.(fields{ii})(:)) > 1
                strptr = 'cstring';
                A.(fields{ii}) = libpointer(strptr, A.(fields{ii}));
            end
        case 'logical'
            if length(A.(fields{ii})(:)) > 1
                strptr = 'int8Ptr';
                A.(fields{ii}) = libpointer(strptr, A.(fields{ii}));
            end
        otherwise
            if length(A.(fields{ii})(:)) > 1
                strptr = [fclass 'Ptr'];
                A.(fields{ii}) = libpointer(strptr, A.(fields{ii}));
            end
    end
end

end