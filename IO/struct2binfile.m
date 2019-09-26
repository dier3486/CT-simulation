function struct2binfile(S, filename, cfgtype)

if nargin<3
    cfgtype = 'xml';
end

cfg = structbincfg(S);
[data, ~] = packstruct(S, cfg);

% output data file
fid = fopen(filename,'w');
fwrite(fid, data, 'uint8');
fclose(fid);

root = struct();
root.(inputname(1)) = cfg;
[pathstr, name, ~] = fileparts(filename);
switch cfgtype
    case 'xml'
        xmlfile = fullfile(pathstr,[name '.xml']);
        struct2xml(root, xmlfile);
    case 'json'
        jsonfile = fullfile(pathstr,[name '.json']);
        jsonwrite(root, jsonfile);
    otherwise
        1;
end