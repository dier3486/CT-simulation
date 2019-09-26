function bincfg = clearbincfg(S, bincfg)
% bincfg = clearbincfg(S, bincfg)

% initial offset
offset_ini = bincfg.offset;
% current offset
offset_cur = offset_ini;

cfgfields = fieldnames(bincfg);
for ifield = 1:length(cfgfields)
    field_ii = cfgfields{ifield};
    switch field_ii
        case 'offset'
            if ischar(bincfg.(field_ii))
                bincfg.(field_ii) = evalcfg(S, bincfg, field);
            end
            if isnan(bincfg.(field_ii))
                bincfg.(field_ii) = offset_cur;
            end
        case {'size', 'number'}
            if ischar(bincfg.(field_ii))
                bincfg.(field_ii) = evalcfg(S, bincfg, field);
            end
        otherwise
            if isstruct(bincfg.(field_ii))
                % recurse
                bincfg.(field_ii) = clearbincfg(S, bincfg.(field_ii));
                offset_cur = offset_cur + bincfg.(field_ii).size * bincfg.(field_ii).number;
            end
    end
end

end

function bincfg = evalcfg(S, bincfg, field)
    statement = bincfg.(field);
    repS = strfind(statement, '$');
    statement(repS) = 'S';
    bincfg.(field) = eval(statement);
end