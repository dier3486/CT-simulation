function [S, bincfg] = sparsepack(data, bincfg, skip)
% sparse bin data to struct
% S = sparsepack(data, bincfg)

if nargin < 3
    skip = 0;
end

% reading file or array
switch class(data)
    case 'uint8'
        filereading = false;
    case 'double'
        filereading = true;
%         pforigin = ftell(data);
    otherwise
        error('error in spasing data!');
end

% decode data size
bincfg.size = decodenumber(bincfg.size, []);
bincfg.number = decodenumber(bincfg.number, []);

if filereading
    % sparse from file
    % datasize left in the file
    datasize = fileposorgtoend(data);
    if ~isavail(bincfg.size)
        % size is unknown
        % try to read first view
        orgnumber = bincfg.number;
        bincfg.number = 1;
%         orgoffest = bincfg.offset;
        [S0, bincfg] = recursesparse(struct(), data, bincfg);
%         packagesize = bincfg.size;
        bincfg.number = orgnumber;
        % filelength
        filelength = datasize / bincfg.size;
        if ~isfield(bincfg, 'filelength') || isavail(bincfg.filelength)
            bincfg.filelength = filelength;
        end
%         bincfg.size = packagesize;
        if ~isavail(bincfg.number) || bincfg.number > (filelength - skip)
            bincfg.number = filelength - skip;
        end
        if bincfg.number == 1 && skip == 0
            S = S0;
            return;
        end
    else
        % filelength
        filelength = datasize / bincfg.size;
        if ~isfield(bincfg, 'filelength') || isavail(bincfg.filelength)
            bincfg.filelength = filelength;
        end
        bincfg.number = min(bincfg.number, filelength - skip);
    end
    bincfg.number = max(bincfg.number, 0);
    % skip views
    fseek(data, bincfg.size*skip, 0);
else
    % sparse from array 
    % empty size or number?
    if ~isavail(bincfg.size) && isavail(bincfg.number)
        % size is unknown
        data = reshape(data, [], bincfg.number + skip);
        bincfg.size = size(data, 1);
    elseif isavail(bincfg.size) && ~isavail(bincfg.number)
        % number is unknown
        data = reshape(data, bincfg.size, []);
        bincfg.number = size(data, 2) - skip;
    elseif ~isavail(bincfg.size) && ~isavail(bincfg.number)
        % both is unknown
%         orgnumber = bincfg.number;
        bincfg.number = 1;
        % try to read first view
        [S0, bincfg] = recursesparse(struct(), data, bincfg);
        % now we kown the size
        data = reshape(data, bincfg.size, []);
        bincfg.number = size(data, 2);
        if bincfg.number == 1 && skip == 0
            S = S0;
            return;
        end
    else
        % both is known
        data = reshape(data, bincfg.size, []);
    end
    % filesize
    if ~isfield(bincfg, 'filesize') || isavail(bincfg.filesize)
        bincfg.filesize = size(data, 2);
    end
    % skip
    data = data(:, skip+1:end);
    if bincfg.number > size(data, 2)
        bincfg.number = size(data, 2);
    end
end

if bincfg.number > 0
    % initial S
    S(bincfg.number) = struct();
else
    S = reshape(struct([]), 1, 0);
end
% recurse sparse
[S, ~] = recursesparse(S, data, bincfg);
% else
%     % to read empty
%     S = reshape(struct([]), 1, 0);
%     S = zerosparse(S, bincfg);
%     1;
% end

end


function [S, bincfg] = recursesparse(S, data, bincfg, offsetshift, datasize, datanumber, Sorig, Sprev)
% sub function recurse
if nargin<4
    offsetshift = 0;
end
if nargin<5
    datasize = bincfg.size;
    datanumber = bincfg.number;
end
if nargin<7
    Sorig = S;
    flag_sorig = true;
else
    flag_sorig = false;
end
if nargin<8
    Sprev = Sorig;
end

% reading file or array
switch class(data)
    case 'uint8'
        filereading = false;
    case 'double'
        filereading = true;
        pforigin = ftell(data);
    otherwise
        error('error in spasing data!');
end

cfgfields = fieldnames(bincfg);
[repS, numS] = size(S);
% offset
offsetshift = offsetshift + bincfg.offset;
for ii = 1:repS
    offsetcount = 0;
    for ifield = 1:length(cfgfields)
        field_ii = cfgfields{ifield};
        cfg_ii = bincfg.(field_ii);
        if ~isstruct(cfg_ii)
            continue
        end
        % numbers
%         cfg_ii.size = decodenumber(S(ii, :), cfg_ii.size, Sprev);
        cfg_ii.size = decodenumber(cfg_ii.size, Sorig, Sprev);
        cfg_ii.number = decodenumber(cfg_ii.number, Sorig, Sprev);
        cfg_ii.offset = decodenumber(cfg_ii.offset, Sorig, Sprev);
        fieldsize = cfg_ii.size * cfg_ii.number;
        if ~isavail(cfg_ii.offset)
            cfg_ii.offset = offsetcount;
        else
            offsetcount = cfg_ii.offset;
        end
        % class
        class_ii = lower(cfg_ii.class);
        if classsize(class_ii)==0 && ~strcmp(class_ii, 'struct')
%             class_ii = decodenumber(S(ii, :), class_ii);
            class_ii = decodenumber(class_ii, Sorig, Sprev);
        end
        cfg_ii.class = class_ii;
        % empty data?
        if cfg_ii.number<=0
            [S(ii, :).(field_ii)] = deal(typecast([], class_ii));
            % to returen bincfg
            if ii == 1
                bincfg.(field_ii) = cfg_ii;
                bincfg.(field_ii).class(1) = upper(bincfg.(field_ii).class(1));
            end
            continue
        end
        switch class_ii
            case 'struct'
                % a sub-struct to recurse
                Srec = struct();
                Srec(cfg_ii.number, numS) = struct();
%                 Srec = recrusesparse(Srec, data, cfg_ii, offsetshift, S(ii, :));
                [Srec, cfg_ii] = recursesparse(Srec, data, cfg_ii, offsetshift, datasize, datanumber, Sorig, S);
                % fieldsize
                if ~isavail(fieldsize)
                    fieldsize = cfg_ii.size * cfg_ii.number;
                end
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
                if filereading
                    fseek(data, offsetshift + cfg_ii.offset, 0);
                    % read data
                    if ifield == length(cfgfields)
                        1;
                    end
                    fielddata = fileblockread(data, fieldsize, datasize, datanumber);
                    % go back
                    fseek(data, pforigin, 'bof');
                else
                    fielddata = data(offsetshift + cfg_ii.offset + (1:fieldsize), 1:datanumber);
                end
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
        
        % to returen bincfg
        if ii == 1
            bincfg.(field_ii) = cfg_ii;
            bincfg.(field_ii).class(1) = upper(bincfg.(field_ii).class(1));
        end

        if flag_sorig
            % updata Sorig and Sprev
            Sorig = S;
            if nargin<8
                Sprev = Sorig;
            end
        end
    end
    if ii == 1 && ~isavail(bincfg.size)
        bincfg.size = offsetcount;
    end
    offsetshift = offsetshift + offsetcount;
end

end


% subfunction to read field data in a file
function fielddata = fileblockread(fid, fieldsize, datasize, datanumber)

if ~isavail(datanumber)
    datanumber = 1;
end

% fielddata = zeros(fieldsize, datanumber, 'uint8');
% pforigin = ftell(fid);
% for ii = 1:datanumber
%     fielddata(:, ii) = fread(fid, fieldsize, 'uint8=>uint8');
%     if ii < datanumber
%         fseek(fid, datasize - fieldsize, 0);
%     end
% end

freadconsol = sprintf('%d*uint8=>uint8', fieldsize);
if isavail(datasize)
    freadskip =  datasize - fieldsize;
else
    freadskip = 0;
end
fielddata = fread(fid, [fieldsize, datanumber], freadconsol, freadskip);

end
