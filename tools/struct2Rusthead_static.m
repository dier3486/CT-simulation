function Rhead = struct2Rusthead_static(A, structname, filename)
% trans matlab structure to RUST head declaration
% Chead = struct2Rusthead_static(A, structname);
% or struct2Rusthead_static(A, structname, filename);
% suggest, struct2Rusthead_static(recon, 'structCTrecon', './declarCTrecon.rs');

if nargin<2 || isempty(structname)
    % default structname to C declaration
    structname = 'mystruct';
end

% recurse the fields to Rust declaration
Cdeclar = Crecurse(A(1), structname, [], 1);

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
C(Cadd).txt = sprintf('pub struct %s {\n', structname);

% loop the fields of A
Afields = fieldnames(A);
for ii = 1:length(Afields)
    txtstr0 = ['\tpub ' Afields{ii} ': '];
    fclass = class(A.(Afields{ii}));
    switch fclass
        case 'struct'
            txtclass = ['struct' structname '_' Afields{ii}];
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
            txtclass = fclass;
    end
    flength = length(A.(Afields{ii})(:));
    if flength > 1
        txtstr = [txtstr0 '[' txtclass ';' num2str(flength) '],\n'];
    else
        txtstr = [txtstr0 txtclass ',\n'];
    end

    if isstruct(A.(Afields{ii}))
        % recurse of the sub struct
        C =  Crecurse(A.(Afields{ii})(1), [structname '_' Afields{ii}], C, level+1);
    end

    % add the txtstr after C.txt
    C(Cadd).txt = [C(Cadd).txt sprintf(txtstr)];
end
% 
C(Cadd).txt = [C(Cadd).txt sprintf('}\n')];

end
