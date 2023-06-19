function A = getstructptr(A)

fields = fieldnames(A);
for ii = 1:length(fields)
    fclass = class(A.(fields{ii}));
    switch fclass
        case 'struct'
            A.(fields{ii}) =  getstructptr(A.(fields{ii}));
        case 'lib.pointer'
            A.(fields{ii}) = get(A.(fields{ii})).Value;
        otherwise
            1;
    end
end

end