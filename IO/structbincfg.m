function bincfg = structbincfg(S, bincfg)
% return the bin configure for a structure

if nargin<2
    bincfg = struct();
    bincfg.offset = 0;
    bincfg.class = 'struct';
    bincfg.size = 0;
end

if ~isstruct(S)
    return
end

% number (length)
bincfg.number = length(S);

% initial offset
offset_ini = bincfg.offset;
% current offset
offset_cur = offset_ini;

% fields
sfields = fieldnames(S);
for ifield = 1: length(sfields)
    % field information
    field_ii = sfields{ifield};
    s_ii = S(1).(field_ii);
    class_ii = class(s_ii);
    
    % to configure
    bincfg.(field_ii).offset = offset_cur;
    bincfg.(field_ii).class = class_ii;
    switch class_ii
        case 'struct'
            bincfg.(field_ii).size = 0;
            % recurse
            bincfg.(field_ii) = structbincfg(s_ii, bincfg.(field_ii));
            offset_cur = offset_cur + bincfg.(field_ii).size * bincfg.(field_ii).number;
        case 'cell'
            % cell is not in supports
            bincfg.(field_ii).number = 0;
            bincfg.(field_ii).size = 0;
            offset_cur = offset_cur + 0;
        otherwise
            bincfg.(field_ii).number = length(s_ii(:));
            bincfg.(field_ii).size =  classsize(class_ii);
            offset_cur = offset_cur + bincfg.(field_ii).size * bincfg.(field_ii).number;
    end
end
bincfg.size = offset_cur - offset_ini;
    
return
