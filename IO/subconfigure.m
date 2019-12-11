function cfg = subconfigure(cfg, isrecurse)
% load configure file in a configure struct
% cfg = subconfigure(cfg, isrecurse)

if nargin<2
    isrecurse = false;
end

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
                if isstruct(cfg.(field_ii){icell}) && isrecurse
                    cfg.(field_ii){icell} = cfgrecurse(cfg.(field_ii){icell});
                end
            end
        end
    end
    if isstruct(cfg.(field_ii)) && isrecurse
        cfg.(field_ii) = cfgrecurse(cfg.(field_ii));
    end
end

end


function subfield = cfgrecurse(subfield)
% recurse

% I know subfield is a strcut
Ns = length(subfield(:));
for ii = 1:Ns
    subfield(ii) = subconfigure(subfield(ii), true);
end

end


function cfg_field = loadsubcfgfile(cfg_field)
% load sub configure file
if exist(cfg_field, 'file')
    [~, ~, cfgext] = fileparts(cfg_field);
    if strcmp(cfgext, '.xml') || strcmp(cfgext, '.json')
        try    
            cfg_field = readcfgfile(cfg_field);
        catch
            % forget it
            1;
        end
    end
end
end
