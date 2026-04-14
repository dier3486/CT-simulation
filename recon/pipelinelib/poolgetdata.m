function [varargout] = poolgetdata(pool, pooldata, Index, getfields, outputtype)
% to get data from pool data (to matlab CPU/gpuArray)

% default pooldata is pool.data
if isempty(pooldata)
    pooldata = pool.data;
end

% the getfields is like, e.g. {'rawdata', 'rawhead_something'}
if ~iscell(getfields)
    getfields = {getfields};
end
Nget = min(size(getfields(:), 1), nargout);


% input type
inputtype = [pool.resourcetype ' ' regexp(pool.bufferresource, '^\w+', 'match', 'once')];
switch inputtype
    case 'Matlab CPU'
        % normal matlab
        inputtype = 1;
    case 'Matlab GPU'
        % gpuArray
        inputtype = 2;
    case 'sharedC CPU'
        % lib.pointer
        inputtype = 3;
    case 'sharedC GPU'
        % lib.pointer on GPU
        inputtype = 4;
    otherwise
        % unknow?
        inputtype = 0;
end

% output type is 
if nargin < 6 || isempty(outputtype)
    outputtype = 0;
elseif ~isnumeric(outputtype)
    switch outputtype
        case 'CPU'
            % normal matlab
            outputtype = 1;
        case 'GPU'
            % gpuArray
            outputtype = 2;
        otherwise
            % default
            outputtype = 0;
    end
end

% Index to pool indexes
index_get =  poolindex(pool, Index);

for iget = 1:Nget
    headkey = regexp(getfields{iget}, '^.+(?=_)', 'match', 'once');
    if ~isempty(headkey)
        if ~isfield(pooldata, headkey) 
            continue; 
        end
        fieldkey = regexp(getfields{iget}, ['(?<=^' headkey '_).+'], 'match', 'once');
        if ~isfield(pooldata.(headkey), fieldkey) 
            continue; 
        end
        % rawhead/imagehead is only on matlab CPU now
        if outputtype == 2
            varargout{iget} = gpuArray(pooldata.(headkey).(fieldkey)(:, index_get));
        else
            varargout{iget} = pooldata.(headkey).(fieldkey)(:, index_get);
        end
    else
        fieldkey = getfields{iget};
        if ~isfield(pooldata, fieldkey) 
            continue; 
        end

        if ~isa(pooldata.(fieldkey), 'lib.pointer') && inputtype > 2
            inputtype_iget = inputtype - 2;
        else
            inputtype_iget = inputtype;
        end

        [datainfo, CUdatainfo] = pooldatainfo(pool.poolsize, pooldata, fieldkey, Index);
        
        switch inputtype_iget * 10 + outputtype
            case {11, 22, 10, 20, 0, 1, 2}
                % CPU -> CPU or GPU -> GPU
                varargout{iget} = pooldata.(fieldkey)(:, index_get);
            case 12
                % CPU -> gpuArray
                varargout{iget} = gpuArray(pooldata.(fieldkey)(:, index_get));
            case 21
                % GPU -> CPU
                varargout{iget} = gather(pooldata.(fieldkey)(:, index_get));
            case {31, 30}
                % lib.pointer -> CPU
                varargout{iget} = libpointer2CPU(pooldata.(fieldkey), datainfo, index_get);
            case 32
                % lib.pointer -> gpuArray
                varargout{iget} = gpuArray(libpointer2CPU(pooldata.(fieldkey), datainfo, index_get));
            case {41, 40}
                % lib.pointer(CUDA) -> CPU
                varargout{iget} = CUDApointer2CPU(pooldata.(fieldkey), datainfo, CUdatainfo);
            case 42
                % lib.pointer(CUDA) -> gpuArray
                varargout{iget} = gpuArray(CUDApointer2CPU(pooldata.(fieldkey), datainfo, CUdatainfo));
            otherwise
                % error
        end

        % varargout{iget} = fieldgetdata(pooldata.(fieldkey), index_get, inputtype, outputtype);
    end

end

end


function data = libpointer2CPU(P, datainfo, index_get)

P.reshape(datainfo.datasize * datainfo.databasis, datainfo.poolsize);
if datainfo.databasis==2 && any(strcmp(P.DataType, {'singlePtr', 'doublePtr'}))
    % complex float
    data = complex(P.Value(1:2:end, index_get), P.Value(2:2:end, index_get));
else
    data = P.Value(:, index_get);
end

end

function data = CUDApointer2CPU(Pgpu, datainfo, CUdatainfo)

if ~libisloaded('CUpipeline')
    % error
    return;
end

Nview = datainfo.index(2) - datainfo.index(1) + 1;
if Nview <=0
    data = zeros(datainfo.datasize * datainfo.databasis, Nview, datainfo.dataclass);
    return;
end
P = libpointer([datainfo.dataclass 'Ptr'], zeros(datainfo.datasize * datainfo.databasis, Nview));

calllib('CUpipeline', 'getpooldata', P, Pgpu, CUdatainfo);

if datainfo.databasis==2 && any(strcmp(datainfo.dataclass, {'single', 'double'}))
    % complex float
    data = complex(P.Value(1:2:end, :), P.Value(2:2:end, :));
else
    data = P.Value;
end

end