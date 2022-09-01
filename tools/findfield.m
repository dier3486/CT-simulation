function [val, address] = findfield(A, fieldname)
% to find a field in struct A

if isfield(A, fieldname)
    val = A(1).(fieldname);
    address = fieldname;
else
    address = '';
    val = [];
    fields = fieldnames(A);
    for ii = 1:length(fields)
        if isstruct(A.(fields{ii}))
            [val, address] = findfield(A(1).(fields{ii}), fieldname);
            if ~isempty(address)
                address = [fields{ii} '.' address];
                break;
            end
        end
    end
end

end