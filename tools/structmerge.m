function A = structmerge(A, B, flag_empty)
% merge struct B to A

if nargin<3
    flag_empty = false;
end

% Afields = fieldnames(A);
Bfields = fieldnames(B);

for ii = 1:length(Bfields)
    field_ii = Bfields{ii};
    if ~isfield(A, field_ii)
        A.(field_ii) = B.(field_ii);
    elseif flag_empty && isempty(A.(field_ii))
        A.(field_ii) = B.(field_ii);
    elseif isstruct(A.(field_ii))
        A.(field_ii) = structmerge(A.(field_ii), B.(field_ii), flag_empty);
    end
end