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
            cfg = jsonread(cfgfile);
        case '.mat'
            cfg = load(cfgfile);
        otherwise
            error(['Unknown file type ',  cfgext]);
    end
else
    error(['File not exist, ', cfgfile]);
end