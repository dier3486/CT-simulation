function S = sparsepack(data, bincfg)
% sparse bin data to struct
% S = sparsepack(data, bincfg)

% reshape data
data = reshape(data, bincfg.size, bincfg.number);
% initial S
S(bincfg.number) = struct();
% I know S is an 1xcfg.number empty struct

% recurse sparse
S = recrusesparse(S, data, bincfg);

end

function S = recrusesparse(S, data, bincfg, offsetshift)

if nargin<4
    offsetshift = 0;
end

cfgfields = fieldnames(bincfg);
[repS, numS] = size(S);
for ii = 1:repS
    for ifield = 1:length(cfgfields)
        field_ii = cfgfields{ifield};
        cfg_ii = bincfg.(field_ii);
        if ~isstruct(cfg_ii)
            continue
        end
        switch cfg_ii.class
            case 'struct'
                % a sub-struct to recurse
                Srec = struct();
                Srec(cfg_ii.number, numS) = struct();
                Srec = recrusesparse(Srec, data, cfg_ii, offsetshift);
                % multi-assigning (row first) with cell
                % cellSrec = mat2cell(Srec, cfg_ii.number, ones(numS,1));
                % multi-assigning (column first) with cell
                cellSrec = mat2cell(Srec.', ones(numS,1), cfg_ii.number);
                [S(ii, :).(field_ii)] = cellSrec{:};
            case 'cell'
                % TBC
                1;
            otherwise
                % sparse a variable (or array) in class 'double', 'single',
                % 'int16' or and so on. 
                fieldsize = cfg_ii.size * cfg_ii.number;
                fielddata = data(offsetshift+cfg_ii.offset+(1:fieldsize), :);
                fielddata = uint8cast(fielddata(:), cfg_ii.class);
                fielddata = reshape(fielddata, cfg_ii.number, []);
                % multi-assigning (row first) with cell
                % celldata = mat2cell(fielddata, cfg_ii.number, ones(numS,1));
                % multi-assigning (column first) with cell
                celldata = mat2cell(fielddata.', ones(numS,1), cfg_ii.number);
                [S(ii, :).(field_ii)] = celldata{:};
                % I know what I did is [a, b, c] = x{:} where x={a, b, c}
        end
    end
    offsetshift = offsetshift + bincfg.size;
end

end



