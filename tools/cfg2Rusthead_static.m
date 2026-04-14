function Rhead = cfg2Rusthead_static(cfg, structname, filename)
% trans matlab structure to RUST head declaration
% Chead = cfg2Rusthead_static(cfg, structname);
% or cfg2Rusthead_static(cfg, structname, filename);
% suggest, cfg2Rusthead_static(cfg, 'rawdata', './declar_rawdata.rs');

if nargin<2 || isempty(structname)
    % default structname to C declaration
    structname = 'mystruct';
end

if ~isstruct(cfg)
    if ~isfile(cfg)
        error('Can not load the cfg.');
    end
    [~, ~, cfgext] = fileparts(cfg);
    switch cfgext
        case '.xml'
            root = myxml2struct(cfg);
        case '.json'
            root = jsonread(cfg);
        otherwise
            error('Can not load the cfg.');
    end
    name = fieldnames(root);
    cfg = root.(name{1});
    structname = name{1};
end

% recurse the fields to Rust declaration
Cdeclar = Crecurse(cfg, structname, [], 1);

% descend sort by .level
[~, levelsort] = sort([Cdeclar(:).level], 'descend');
Cdeclar = Cdeclar(levelsort);

% return the txt
Rhead = [Cdeclar(:).txt];

% write to filename
if exist('filename', 'var')
    fid = fopen(filename, 'w');
    fwrite(fid, Rhead);
    fclose(fid);
end

end


function C = Crecurse(A, structname, C, level)

% add new element after the end of C
Cadd = length(C(:)) + 1;

% level
C(Cadd).level = level;
C(Cadd).structname = structname;  % useless

% kick off txt
C(Cadd).txt = '';
% repr line
C(Cadd).txt = [C(Cadd).txt sprintf('#[repr(C, packed(1))]\n')];
% derive line
C(Cadd).txt = [C(Cadd).txt sprintf('#[derive(Default, Clone, Copy, Debug)]\n')];
C(Cadd).txt = [C(Cadd).txt sprintf('pub struct %s {\n', structname)];

offset_cum = 0;
resv = 0;
% loop the fields of A
Afields = fieldnames(A);
for ii = 1:length(Afields)
    if ~isstruct(A.(Afields{ii}))
        continue;
    end
    fclass = A.(Afields{ii}).class;
    offset = A.(Afields{ii}).offset;
    number = A.(Afields{ii}).number;
    size = A.(Afields{ii}).size;
    if isnumeric(offset)
        if isnan(offset_cum)
            offset_cum = offset;
        elseif offset > offset_cum
            % fill a reserve
            txtresv = ['\t reserve_' num2str(resv) ': [u8;' num2str(offset - offset_cum) '],\n'];
            C(Cadd).txt = [C(Cadd).txt sprintf(txtresv)];
            resv = resv + 1;
            offset_cum = offset;
        end
    end
    if isnumeric(number) && isnumeric(size)
        offset_cum = offset_cum + number*size;
    else
        offset_cum = nan;
    end
    
    torecurse = false;
    switch lower(fclass)
        case 'struct'
            txtclass = ['struct' structname '_' Afields{ii}];
            torecurse = true;
        case 'char'
            txtclass = 'char';
        case 'uint8'
            txtclass = 'u8';
        case 'uint16'
            txtclass = 'u16';
        case 'uint32'
            txtclass = 'u32';
        case 'uint64'
            txtclass = 'u64';
        case 'int8'
            txtclass = 'i8';
        case 'int16'
            txtclass = 'i16';
        case 'int32'
            txtclass = 'i32';
        case 'int64'
            txtclass = 'i64';
        case 'single'
            txtclass = 'f32';
        case 'double'
            txtclass = 'f64';
            % no comlex plz
        case 'logical'
            txtclass = 'bool';
        otherwise
            txtclass = 'NA';
    end

    txtstr0 = ['\t' Afields{ii} ': '];
    if ~isnumeric(number)
        txtstr = [txtstr0 '[' txtclass ';NA],\n'];
    elseif number > 1
        txtstr = [txtstr0 '[' txtclass ';' num2str(number) '],\n'];
    else
        txtstr = [txtstr0 txtclass ',\n'];
    end

    if torecurse
        % recurse of the sub struct
        C =  Crecurse(A.(Afields{ii}), [structname '_' Afields{ii}], C, level+1);
    end

    % add the txtstr after C.txt
    C(Cadd).txt = [C(Cadd).txt sprintf(txtstr)];
end
% 
C(Cadd).txt = [C(Cadd).txt sprintf('}\n')];

end
