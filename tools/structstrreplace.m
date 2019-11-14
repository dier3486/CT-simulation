function S = structstrreplace(S, torep, repby)
% replace string torep by repby

Sclass = class(S);
sizeS = size(S(:), 1);
switch Sclass
    case 'char'
        % replace string
        S = strrep(S, torep, repby);
    case 'cell'
        % loop cell to recurse
        for ii = 1:sizeS
            S{ii} = structstrreplace(S{ii}, torep, repby);
        end
    case 'struct'
        % loop struct to recurse
        if sizeS>1
            for ii = 1:sizeS
                S(ii) = structstrreplace(S(ii), torep, repby);
            end
        else
            fields = fieldnames(S);
            for ifield = 1:length(fields)
                S.(fields{ifield}) = structstrreplace(S.(fields{ifield}), torep, repby);
            end
        end
    otherwise
        % do nothing
        1;
end