function cfg = subconfigure(cfg)
% load configure file in a configure struct

fields = fieldnames(cfg);
for ii = 1:length(fields)
    field_ii = fields{ii};
    if ischar(cfg.(field_ii))
        % a string
        cfg.(field_ii) = loadsubcfgfile(cfg.(field_ii));
    elseif iscell(cfg.(field_ii))
        % a cell
        for icell = 1:length(cfg.(field_ii)(:))
            % loop the cell
            if ischar(cfg.(field_ii){icell})
                cfg.(field_ii){icell} = loadsubcfgfile(cfg.(field_ii){icell});
            end
        end
    end
end
end


function cfg_field = loadsubcfgfile(cfg_field)
% load sub configure file
if exist(cfg_field, 'file')
    try
        cfg_field = readcfgfile(cfg_field);
    catch
        % forget it
        1;
    end
end
end
