function [r, c] = structcompare(A, B, flag_size, c)
% to return A==B,
%   r = structcompare(A, B);
% it is a debug tool, no GPU

if nargin < 3
    flag_size = false;
end
if nargin < 4
    c = 'root';
end

Afields = fieldnames(A);
Bfields = fieldnames(B);
if length(Afields) ~= length(Bfields)
    r = false;
    return;
end

Asize = size(A);
Bsize = size(B);
if flag_size
    try
        s = ~all(Asize==Bsize);
    catch
        r = false;
        return
    end
    if s
        r = false;
        return
    elseif prod(Asize)==0
        r = true;
        return
    end
else
    if prod(Asize) ~= prod(Bsize)
        r = false;
        return
    elseif prod(Asize)==0
        r = true;
        return
    end
end

r = true;
for ii = 1:length(Afields)
    field_ii = Afields{ii};
    if ~isfield(B, field_ii)
        r = false;
        c = [c '.' field_ii];
        return;
    end
    classA = class(A(1).(field_ii));
    classB = class(B(1).(field_ii));
    if ~strcmp(classA, classB)
        r = false;
        c = [c '.' field_ii];
        return;
    end
    for jj = 1 : prod(Asize)
        if prod(Asize) > 1
            strjj = sprintf('(%d)', jj);
        else
            strjj = '';
        end
        switch classA
            case 'struct'
                % to recurese
                [r, cr] = structcompare(A(jj).(field_ii), B(jj).(field_ii), flag_size, [c strjj '.' field_ii]);
                if ~r
                    c = cr;
                    return
                end
            case 'cell'
                [r, cr] = cellrecurse(A(jj).(field_ii), B(jj).(field_ii), flag_size, [c strjj '.' field_ii]);
                if ~r
                    c = cr;
                    return
                end
                %error('Sorry, no cell plz.');
            otherwise
                r = comparevalue(A(jj).(field_ii), B(jj).(field_ii), flag_size);
                if ~r
                    c = [c strjj '.' field_ii];
                    return
                end
        end

    end
end

if nargin < 4 && r
    c = '';
end

end


function r = comparevalue(A_field, B_field, flag_size)

Afsize = size(A_field);
Bfsize = size(B_field);
if flag_size
    try
        s = ~all(Afsize==Bfsize);
    catch
        r = false;
        return
    end
    if s
        r = false;
        return
    elseif prod(Afsize)==0
        r = true;
        return
    end
else
    if prod(Afsize) ~= prod(Bfsize)
        r = false;
        return
    elseif prod(Afsize)==0
        r = true;
        return
    end
end
r = all(isnan(A_field(:)) == isnan(B_field(:)));
fillconst = cast(0, class(A_field));
r = r & all(fillmissing(A_field(:), 'constant', fillconst) == ...
    fillmissing(B_field(:), 'constant', fillconst));

end

function [r, c] = cellrecurse(A, B, flag_size, c)

Asize = size(A);
Bsize = size(B);
if flag_size
    try
        s = ~all(Asize==Bsize);
    catch
        r = false;
        return
    end
    if s
        r = false;
        return
    elseif prod(Asize)==0
        r = true;
        return
    end
else
    if prod(Asize) ~= prod(Bsize)
        r = false;
        return
    elseif prod(Asize)==0
        r = true;
        return
    end
end


r = true;
for jj = 1 : prod(Asize)
    strjj = sprintf('{%d}', jj);
    classA = class(A{jj});
    classB = class(B{jj});
    if ~strcmp(classA, classB)
        r = false;
        c = [c strjj];
        return;
    end
    switch classA
        case 'struct'
            % to recurese
            [r, cr] = structcompare(A{jj}, B{jj}, flag_size, [c strjj]);
            if ~r
                c = cr;
                return
            end
        case 'cell'
            % cell in cell
            [r, cr] = cellrecurse(A{jj}, B{jj}, flag_size, [c strjj]);
            if ~r
                c = cr;
                return
            end
            %error('Sorry, no cell plz.');
        otherwise
            r = comparevalue(A{jj}, B{jj}, flag_size);
            if ~r
                c = [c strjj];
                return
            end
    end

end

end
