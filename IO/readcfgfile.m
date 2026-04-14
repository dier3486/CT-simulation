function cfg = readcfgfile(cfgfile)

cfgfile = char(cfgfile);
if exist(cfgfile, 'file')
    [~, ~, cfgext] = fileparts(cfgfile);
    switch cfgext
        case '.xml'
            root = myxml2struct(cfgfile);
            name = fieldnames(root);
            cfg = root.(name{1});
        case '.json'
            root = jsonread(cfgfile);
            name = fieldnames(root);
            cfg = root.(name{1});
            % decode 'nan', 'inf', '-inf'
            cfg = checkcharnum(cfg);
        case '.mat'
            cfg = load(cfgfile);
        otherwise
            error(['Unknown file type ',  cfgext]);
    end
else
    error(['File not exist, ', cfgfile]);
end

end


function A = checkcharnum(A, specialnum)

if nargin<2
    specialnum = {'nan', 'inf', '-inf'};
elseif ~iscell(specialnum)
    specialnum = {specialnum};
end

classA = class(A);
switch classA
    case 'struct'
        if length(A)>1
            for ii = 1:length(A(:))
                A(ii) = checkcharnum(A(ii), specialnum);
            end
        else
            fields = fieldnames(A);
            for ifield = 1:length(fields)
                A.(fields{ifield}) = checkcharnum(A.(fields{ifield}), specialnum);
            end
        end
    case 'cell'
        for icell = 1:length(A(:))
            A{icell} = checkcharnum(A{icell}, specialnum);
        end
    case {'char', 'string'}
        if any(strcmpi(specialnum, A))
            A = eval(A);
        end
    otherwise
        0;
end

end