function [nextdata, writenum2] = pooldatacopy(currpool, currdata, nextpool, nextdata, writenum, ...
    copyfields, force_flag, method)
% to copy the data from currunt data to next data in pipeline, circulated patch
%   nextdata = pooldatacopy(currpool, currdata, nextpool, nextdata, writenum, copyfields, force_flag);
% or
%   nextdata = pooldatacopy(currpool, currdata, nextpool, nextdata, writenum);
% It will be moved to pooldatacopy as an inner subfucntion

if nargin < 7
    force_flag = false;
end

if nargin < 6 || isempty(copyfields)
    if force_flag
        copyfields = currpool.datafields;
    else
        copyfields = nextpool.datafields;
    end
end

if nargin < 8
    method = 'overwrite';
end

[ReadPoint, WritePoint, writenum2] = checkwritenumber(currpool, nextpool, currpool.ReadPoint, nextpool.WritePoint, writenum);
if writenum2 ~= writenum
    % play again
    [ReadPoint, WritePoint, writenum2] = checkwritenumber(currpool, nextpool, ReadPoint, WritePoint, writenum2);
end
if writenum2 <= 0
    writenum2 = 0;
    return;
end

nextdata = poolhardcopy2(nextdata, currdata, nextpool, currpool, WritePoint, ReadPoint, ...
    writenum2, copyfields, method);
end


function [ReadPoint, WritePoint, writenum] = checkwritenumber(currpool, nextpool, ReadPoint, WritePoint, writenum)

if currpool.circulatemode
    ReadPoint = mod(ReadPoint - 1, currpool.poolsize) + 1;
else
    if ReadPoint < 1
        % skip the negative indeces
        writenum = writenum - ReadPoint - 1;
        WritePoint = WritePoint + ReadPoint + 1;
        ReadPoint = 1;
    end
    if ReadPoint + writenum - 1 > currpool.poolsize
        % skip the outof range indeces
        writenum = currpool.poolsize - ReadPoint + 1;
    end
end
if nextpool.circulatemode
    WritePoint = mod(WritePoint - 1, nextpool.poolsize) + 1;
else
    if WritePoint < 1
        % skip the negative indeces
        writenum = writenum - WritePoint - 1;
        ReadPoint = ReadPoint + WritePoint + 1;
        WritePoint = 1;
    end
    if WritePoint + writenum - 1 > nextpool.poolsize
        % skip the outof range indeces
        writenum = nextpool.poolsize - WritePoint + 1;
    end
end

end


function outdata = poolhardcopy2(outdata, indata, outpool, inpool, PointOut, PointIn, copynumber, copyfields, method)

if isempty(copyfields)
    copyfields = fieldnames(indata);
end
if isempty(method)
    method = 'overwrite';
    % or 'cum'
end

% to recurse the struct
for ii = 1:length(copyfields)
    if isempty(copyfields{ii}) || ~isfield(indata, copyfields{ii})
        % pass the {''} and the not exist fields
        continue;
    end
    if size(indata.(copyfields{ii}), 2) == 1 && isstruct(indata.(copyfields{ii}))
        % recurse the struct
        if ~isfield(outdata, copyfields{ii})
            outdata.(copyfields{ii}) = struct();
        end
        outdata.(copyfields{ii}) = ...
            poolhardcopy2(outdata.(copyfields{ii}), indata.(copyfields{ii}), outpool, inpool, PointOut, PointIn, ...
            copynumber, [], method);
    end
end

% index
Index_in = [PointIn, PointIn + copynumber - 1];
Index_out = [PointOut, PointOut + copynumber - 1];
index_in = poolindex(inpool, Index_in);
index_out = poolindex(outpool, Index_out);

% loop the fields
for ii = 1:length(copyfields)
    % pass the {''} and the not exist fields
    if isempty(copyfields{ii}) || ~isfield(indata, copyfields{ii})
        continue;
    end
    % pass the structs (which were recursed)
    if size(indata.(copyfields{ii}), 2) == 1 && isstruct(indata.(copyfields{ii}))
        % the struct was recursed
        continue;
    end
    % check the empty or missed fields in outdata
    if ~isfield(outdata, copyfields{ii}) || all(size(outdata.(copyfields{ii})) == 0)
        % shall not happen on deployment
        outdata.(copyfields{ii}) = zeros(size(indata.(copyfields{ii}), 1), 0, 'like', indata.(copyfields{ii}));
        outdata.(copyfields{ii})(:, index_out) = indata.(copyfields{ii})(:, index_in);
        continue;
    end
    % the lib.pointers
    isOutPtr = isa(outdata.(copyfields{ii}), 'lib.pointer');
    isInPtr = isa(indata.(copyfields{ii}), 'lib.pointer');
    % check the lacking size of the outdata
    if size(outdata.(copyfields{ii}), 2) < max(index_out) && ~isOutPtr
        outdata.(copyfields{ii})(:, max(index_out)) = 0;
    end
    % while the outdata isreal only copy the real part
    if isfield(indata, [copyfields{ii}, '_databasis'])
        inreal = indata.([copyfields{ii}, '_databasis']) ~= 2;
    elseif ~isInPtr
        inreal = isreal(indata.(copyfields{ii}));
    else
        inreal = true;
    end
    if isfield(outdata, [copyfields{ii}, '_databasis'])
        outreal = outdata.([copyfields{ii}, '_databasis']) ~= 2;
    elseif ~isOutPtr
        outreal = isreal(outdata.(copyfields{ii}));
    else
        outreal = inreal;
    end
    if ~isOutPtr
        dataget = poolgetdata(inpool, indata, Index_in, copyfields{ii}, 0);
        % commit the copy operator
        switch lower(method)
            case 'overwrite'
                if outreal && ~inreal
                    outdata.(copyfields{ii})(:, index_out) = real(dataget);
                % elseif ~outreal && inreal
                else
                    outdata.(copyfields{ii})(:, index_out) = dataget;
                end
            case 'cum'
                if outreal && ~inreal
                    outdata.(copyfields{ii})(:, index_out) = ...
                        outdata.(copyfields{ii})(:, index_out) + real(dataget);
                else
                    outdata.(copyfields{ii})(:, index_out) = ...
                        outdata.(copyfields{ii})(:, index_out) + dataget;
                end
            case 'clear'
                outdata.(copyfields{ii})(:, index_out) = 0;
            otherwise
                % do nothing
                0;
                % We can support function handles here, later.
        end
        % to fix the auto complex -> real bug (need not in deployment codes)
        if ~outreal && isreal(outdata.(copyfields{ii}))
            % it is a bug of matlab, I can not even avoid it by using cast or more if-else.
            outdata.(copyfields{ii}) = complex(outdata.(copyfields{ii}));
        end
        % I know the GPU/CPU can be autoly gathered or put.
    elseif ~isInPtr
        outdata = poolputdata(outpool, outdata, Index_out, indata.(copyfields{ii}), copyfields{ii});
    else  % Ptr to Ptr
        % copy prepare
        [~, CUindatainfo] = pooldatainfo(inpool.poolsize, indata, copyfields{ii}, Index_in);
        [~, CUoutdatainfo] = pooldatainfo(outpool.poolsize, outdata, copyfields{ii}, Index_out);
        % copy kind
        isOutGPU = strcmpi(outpool.bufferresource(1:3), 'GPU');
        isInGPU = strcmpi(inpool.bufferresource(1:3), 'GPU');
        % copyKind = 0 : Host to Host, 1: Host to Device, 2: Device to Host, 3: Device to Device
        copyKind = isInGPU * 2 + isOutGPU;
        switch lower(method)
            case 'overwrite'
                calllib('CUpipeline', 'pooldatacopy', ...
                    oudata.(copyfields{ii}), indata.(copyfields{ii}), CUoutdatainfo, CUindatainfo, copyKind);
            case 'cum'
                if copyKind == 0
                    outdata.(copyfields{ii}).Value(:, index_out) = outdata.(copyfields{ii}).Value(:, index_out) + ...
                        indata.(copyfields{ii}).Value(:, index_in); 
                else
                    error('Can not employ any extra operator in CUDA memory copy!');
                end
            case 'clear'
                if isOutGPU
                    calllib('CUpipeline', 'pooldataclear', outdata.(copyfields{ii}), CUdatainfo);
                else
                    outdata.(copyfields{ii}).Value(:, index_out) = 0;
                end
            otherwise
                % do nothing
                0;
        end
    end
end

end
