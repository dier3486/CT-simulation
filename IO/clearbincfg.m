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
%             if ischar(bincfg.(field_ii))
%                 bincfg = evalcfg(S, bincfg, field_ii);
%             end
            bincfg.(field_ii) = decodenumber(S, bincfg.(field_ii));
            if isnan(bincfg.(field_ii))
                bincfg.(field_ii) = offset_cur;
            end
        case {'size', 'number'}
%             if ischar(bincfg.(field_ii))
%                 bincfg = evalcfg(S, bincfg, field_ii);
%             end
            bincfg.(field_ii) = decodenumber(S, bincfg.(field_ii));
            if isnan(bincfg.(field_ii))
                bincfg.(field_ii) = [];
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

function r = decodenumber(S, c)
% explain the numers in cfg_ii.size and cfg_ii.number
    if isnumeric(c)
        r = c;
    elseif isempty(c)
        r = [];
    elseif ischar(c)
%         c(c=='$') = 'S';
        c = strrep(c, '$', 'S');
        r = eval(c);
    else
        r = nan;
    end

end

% function bincfg = evalcfg(S, bincfg, field)
%     statement = bincfg.(field);
%     repS = strfind(statement, '$');
%     statement(repS) = 'S';
%     bincfg.(field) = eval(statement);
% end