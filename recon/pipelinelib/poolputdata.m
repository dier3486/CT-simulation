function pooldata = poolputdata(pool, pooldata, Index, dataIn, putfield)
% to put (matlab) data to pool data

% default pooldata is pool.data
if isempty(pooldata)
    pooldata = pool.data;
end

% output type
outputtype = [pool.resourcetype ' ' regexp(pool.bufferresource, '^\w+', 'match', 'once')];

% Index to pool indexes
index_put =  poolindex(pool, Index);

% the putfield is like, e.g. 'rawdata' or 'rawhead_something'.
headkey = regexp(putfield, '^.+(?=_)', 'match', 'once');

if ~isempty(headkey)
    fieldkey = regexp(putfield, ['(?<=^' headkey '_).+'], 'match', 'once');
    % rawhead/imagehead is only on matlab CPU now
    pooldata.(headkey).(fieldkey)(:, index_put) = gather(dataIn);
else
    fieldkey = putfield;
    % copy prepare
    [datainfo, CUdatainfo] = pooldatainfo(pool.poolsize, pooldata, fieldkey, Index);

    switch outputtype
        case 'Matlab CPU'
            % CPU/gpuArray -> CPU
            pooldata.(fieldkey)(:, index_put) = gather(dataIn);
        case 'Matlab GPU'
            % CPU/gpuArray -> gpuArray
            pooldata.(fieldkey)(:, index_put) = gpuArray(dataIn);
        case 'sharedC CPU'
            % CPU/gpuArray -> lib.pointer
            data2libpointer(pooldata.(fieldkey), gather(dataIn), datainfo, index_put);
        case 'sharedC GPU'
            % CPU/gpuArray -> lib.pointer(CUDA)
            data2CUDApointer(pooldata.(fieldkey), gather(dataIn), datainfo, CUdatainfo);
        otherwise
            % CPU -> CPU
            pooldata.(fieldkey)(:, index_put) = dataIn;
    end
end



end


function data2libpointer(P, dataIn, datainfo, index_put)

P.reshape(datainfo.datasize * datainfo.databasis, datainfo.poolsize);
if datainfo.databasis==2 && any(strcmp(P.DataType, {'singlePtr', 'doublePtr'}))
    % complex float
    P.Value(1:2:end, index_put) = real(dataIn);
    P.Value(2:2:end, index_put) = imag(dataIn);
else
    P.Value(:, index_put) = dataIn;
end

end

function data2CUDApointer(P, dataIn, datainfo, CUdatainfo)

if ~libisloaded('CUpipeline')
    % error
    return;
end

Nview = datainfo.index(2) - datainfo.index(1) + 1;
if Nview <=0
    return;
end

if datainfo.databasis==2 && any(strcmp(datainfo.dataclass, {'single', 'double'})) 
    Pcpu = libpointer([datainfo.dataclass 'Ptr'], complex2float(dataIn));
else
    Pcpu = libpointer([datainfo.dataclass 'Ptr'], dataIn);
end

calllib('CUpipeline', 'putpooldata', P, Pcpu, CUdatainfo);

end