function A = structcopy(A, B, flag_recurse)
% copy struct B to A, but skip the fields not existed in A
% A = structcopy(A, B, flag_recurse);
% FLAG:
%   flag_recurse        false:  no recurse
%                       true:   recurse the sub structs (default)

if nargin<3
    flag_recurse = true;
end

if isempty(A)
    A = struct();
end
if isempty(B)
    B = struct();
end

Afields = fieldnames(A);
Bfields = fieldnames(B);

for ii = 1:length(Bfields)
    field_ii = Bfields{ii};
    if ~isincell(Afields, field_ii)
        continue;
    end
    if isstruct(A.(field_ii)) && flag_recurse
        A.(field_ii) = structcopy(A.(field_ii), B.(field_ii), flag_recurse);
    else
        A.(field_ii) = B.(field_ii);
    end
end