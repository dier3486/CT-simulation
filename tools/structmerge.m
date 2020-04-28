function A = structmerge(A, B, flag_empty, flag_recurse)
% merge struct B to A, but not overwriting the fields already existed in A
% A = structmerge(A, B, flag_empty, flag_recurse);
% FLAG:
%   flag_empty          false:  skip the empty fields in A (default)
%                       true:   overwrite the empty fields in A
%   flag_recurse        false:  no recurse
%                       true:   recurse the sub structs (default)

if nargin<3
    flag_empty = false;
end
if nargin<4
    flag_recurse = true;
end

if isempty(A)
    A = struct();
end
if isempty(B)
    B = struct();
end

% Afields = fieldnames(A);
Bfields = fieldnames(B);

for ii = 1:length(Bfields)
    field_ii = Bfields{ii};
    if ~isfield(A, field_ii)
        A.(field_ii) = B.(field_ii);
    elseif flag_empty && isempty(A.(field_ii))
        A.(field_ii) = B.(field_ii);
    elseif isstruct(A.(field_ii)) && flag_recurse
        A.(field_ii) = structmerge(A.(field_ii), B.(field_ii), flag_empty, flag_recurse);
    end
end