function bincfg = clearbincfg(bincfg, Sorig, Sprev, Scurr)
% subfucntion only used in packstruct.m
% bincfg = clearbincfg(bincfg, S)

if nargin < 3
    Sprev = Sorig;
end
if nargin < 4
    Scurr = Sprev;
end

% initial offset
offset_ini = bincfg.offset;
% current offset
offset_cur = offset_ini;

cfgfields = fieldnames(bincfg);
for ifield = 1:length(cfgfields)
    field_ii = cfgfields{ifield};
    switch field_ii
        case 'offset'
            bincfg.(field_ii) = decodenumber(bincfg.(field_ii), Sorig, Sprev);
            if size(bincfg.(field_ii)(:), 1) > 1
                bincfg.(field_ii) = bincfg.(field_ii)(1);
            end
            if isnan(bincfg.(field_ii))
                bincfg.(field_ii) = offset_cur;
            end
        case {'size', 'number'}
            bincfg.(field_ii) = decodenumber(bincfg.(field_ii), Sorig, Sprev);
            if size(bincfg.(field_ii)(:), 1) > 1
                bincfg.(field_ii) = bincfg.(field_ii)(1);
            end
            if isnan(bincfg.(field_ii))
                bincfg.(field_ii) = [];
            end
        case 'class'
            % do nothing
            % class can not be flexible
            1;
        otherwise
            if isstruct(bincfg.(field_ii))
                % recurse
                if isfield(Scurr, field_ii) && ~isempty(Scurr.(field_ii))
                    bincfg.(field_ii) = clearbincfg(bincfg.(field_ii), Sorig, Scurr, Scurr.(field_ii)(1));
                else
                    bincfg.(field_ii) = clearbincfg(bincfg.(field_ii), Sorig, Scurr, struct());
                end
                offset_cur = offset_cur + double(bincfg.(field_ii).size) * double(bincfg.(field_ii).number);
            end
    end
end

if ~isavail(bincfg.size)
    bincfg.size = offset_cur;
end
if ~isavail(bincfg.number)
    bincfg.number = length(Scurr(:));
end

end
