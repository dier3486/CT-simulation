function S = sparsepack2(data, bincfg)
% sparse bin data to struct
% S = sparsepack2(data, bincfg)

% reshape data
bincfg.size = decodenumber([], bincfg.size);
bincfg.number = decodenumber([], bincfg.number);

% empty size or number?
if isempty(bincfg.size) && ~isempty(bincfg.number)
    % size is unknown
    data = reshape(data, bincfg.size, bincfg.number);
    bincfg.size = size(data, 1);
elseif isempty(bincfg.number) && ~isempty(bincfg.size)
    % number is unknown
    data = reshape(data, bincfg.size, bincfg.number);
    bincfg.number = size(data, 2);
else
    % both is unknown
    data = data(:);
    [bincfg.size, bincfg.number] = size(data);
    % try to read first view
    [~, count1] = recursesparse(struct(), data, bincfg);
    bincfg.size = count1;
    data = reshape(data, bincfg.size, []);
    bincfg.number = size(data, 2);
end

% initial S
S(bincfg.number) = struct();
% I know S is an 1xcfg.number empty struct

% recurse sparse
S = recursesparse(S, data, bincfg);

end


function [S, offsetcount] = recursesparse(S, data, bincfg, offsetshift, Sprev)
% sub function recurse
if nargin<4
    offsetshift = 0;
end
if nagrin<5
    Sprev = S;
end

cfgfields = fieldnames(bincfg);
[repS, numS] = size(S);
for ii = 1:repS
    offsetcount = bincfg.offset;
    for ifield = 1:length(cfgfields)
        field_ii = cfgfields{ifield};
        cfg_ii = bincfg.(field_ii);
        if ~isstruct(cfg_ii)
            continue
        end
        % numbers
        cfg_ii.size = decodenumber(S(ii, :), cfg_ii.size, Sprev);
        cfg_ii.number = decodenumber(S(ii, :), cfg_ii.number, Sprev);
        cfg_ii.offset = decodenumber(S(ii, :), cfg_ii.offset, Sprev);
        fieldsize = cfg_ii.size * cfg_ii.number;
        if isnan(cfg_ii.offset)
            cfg_ii.offset = offsetcount;
        else
            offsetcount = cfg_ii.offset;
        end
        % class
        class_ii = lower(cfg_ii.class);
        if classsize(class_ii)==0 && ~strcmp(class_ii, 'struct')
            class_ii = decodenumber(S(ii, :), class_ii);
        end
        if cfg_ii.number<=0
            [S(ii, :).(field_ii)] = deal(typecast([], class_ii));
            continue
        end
        switch class_ii
            case 'struct'
                % a sub-struct to recurse
                Srec = struct();
                Srec(cfg_ii.number, numS) = struct();
                Srec = recrusesparse(Srec, data, cfg_ii, offsetshift, S(ii, :));
                % multi-assigning (row first) with cell
                cellSrec = mat2cell(Srec, cfg_ii.number, ones(numS,1));
                % multi-assigning (column first) with cell
                % cellSrec = mat2cell(Srec.', ones(numS,1), cfg_ii.number);
                [S(ii, :).(field_ii)] = cellSrec{:};
            case 'cell'
                % TBC
                1;
            otherwise
                % sparse a variable (or array) in class 'double', 'single',
                % 'int16' or and so on.
                % get data
                fielddata = data(offsetshift + cfg_ii.offset + (1:fieldsize), :);
                fielddata = uint8cast(fielddata(:), class_ii);
                fielddata = reshape(fielddata, cfg_ii.number, []);
                % multi-assigning (row first) with cell
                celldata = mat2cell(fielddata, cfg_ii.number, ones(numS,1));
                % multi-assigning (column first) with cell
                % celldata = mat2cell(fielddata.', ones(numS,1), cfg_ii.number);
                [S(ii, :).(field_ii)] = celldata{:};
                % I know what I did is [a, b, c] = x{:} where x={a, b, c}
        end
        offsetcount = offsetcount + fieldsize;
    end
    offsetshift = offsetshift + bincfg.size;
end

end


function r = decodenumber(S, c, Sprev)
% explain the numers in cfg_ii.size and cfg_ii.number
    if isnumeric(c)
        r = c;
    elseif isempty(c)
        r = [];
    elseif ischar(c)
%         c(c=='$') = 'S';
        c = regexprep(c, '\$\$' , 'Sprev');
        c = regexprep(c, '\$' , 'S');
        try
            r = eval(c);
        catch
            r = c;
        end
    else
        r = nan;
    end
    % to double
    if isnumeric(r)
        r = double(r);
    end

end
