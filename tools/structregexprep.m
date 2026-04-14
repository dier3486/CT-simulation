function S = structregexprep(S, torep, repby, root)
% replace the text in S matched with expression torep by repby

if nargin<4
    root = S;
end

Sclass = class(S);
sizeS = size(S(:), 1);
switch Sclass
    case 'char'
        % replace string
        try
            S = regexprep(S, torep, repby);
        catch
            % normaly the field defined repby not exist in root
            S = regexprep(S, torep, '');
            if any(strcmpi(S, {'/', '\'}))
                % do not let S be '/'
                S = '';
            end
        end
    case 'cell'
        % loop cell to recurse
        for ii = 1:sizeS
            S{ii} = structregexprep(S{ii}, torep, repby, root);
        end
    case 'struct'
        % loop struct to recurse
        if sizeS>1
            for ii = 1:sizeS
                S(ii) = structregexprep(S(ii), torep, repby, root);
            end
        else
            fields = fieldnames(S);
            for ifield = 1:length(fields)
                S.(fields{ifield}) = structregexprep(S.(fields{ifield}), torep, repby, root);
            end
        end
    otherwise
        % do nothing
        1;
end