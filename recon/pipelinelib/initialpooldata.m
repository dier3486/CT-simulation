function pooldata = initialpooldata(pooldata, objecttype, poolsize, datasize, databasis, gpuDevice, resourcetype)
% to initial the pool.data

if isempty(pooldata)
    pooldata = struct();
end

% head
switch objecttype
    case 'rawdata'
        pooldata.rawhead = struct();
    case {'image', 'sparseimage'}
        pooldata.imagehead = struct();
    otherwise
        pooldata.([objecttype 'head']) = struct();
end

% data
switch objecttype
    case {'rawdata', 'image'}
        % datasize and databasis method
        pooldata.([objecttype '_datasize']) = datasize;
        pooldata.([objecttype '_databasis']) = databasis;
        % initial the pool data
        switch resourcetype
            case 'Matlab'
                if gpuDevice > 0
                    if databasis==2
                        % complex
                        pooldata.(objecttype) = complex(zeros(datasize, poolsize, 'single', 'gpuArray'));
                    else % real or databasis>2
                        pooldata.(objecttype) = zeros(datasize*databasis, poolsize, 'single', 'gpuArray');
                    end
                else
                    if databasis==2
                        % complex
                        pooldata.(objecttype) = complex(zeros(datasize, poolsize, 'single'));
                    else % real or databasis>2
                        pooldata.(objecttype) = zeros(datasize*databasis, poolsize, 'single');
                    end
                end
            case 'sharedC'
                if gpuDevice > 0
                    % CUDA
                    datapoolsize = uint64(datasize * databasis * poolsize * 4);
                    if isfield(pooldata, objecttype) && isa(pooldata.(objecttype), 'lib.pointer')
                        % free it before renew
                        calllib('CUpipeline', 'cudaworkspaceFree', pooldata.(objecttype));
                    end
                    pooldata.(objecttype) = calllib('CUpipeline', 'cudaworkspaceNew', datapoolsize);
                else
                    % lib.pointer
                    pooldata.(objecttype) = libpointer('singlePtr', zeros(datasize*databasis, poolsize));
                end
            otherwise
        end
    case 'sparseimage'
        % datasize and databasis method
        pooldata.image_datasize = datasize;
        pooldata.image_databasis = databasis;
        % pooldata.image_datannz = 0;
        % initial the pool data
        switch resourcetype
            case 'Matlab'
                if gpuDevice > 0
                    if databasis==2
                        % complex
                        pooldata.image = complex(zeros(datasize, poolsize, 'single', 'gpuArray'));
                    else % rela or databasis>2
                        pooldata.image = zeros(datasize*databasis, poolsize, 'single', 'gpuArray');
                    end
                    pooldata.Sreconact = false(datasize, poolsize, 'gpuArray');
                    % or to use CSR
                    pooldata.CSR = zeros(datasize, poolsize, 'int32', 'gpuArray');
                else
                    if databasis==2
                        % complex
                        pooldata.image = complex(zeros(datasize, poolsize, 'single'));
                    else % real or databasis>2
                        pooldata.image = zeros(datasize*databasis, poolsize, 'single');
                    end
                    % Sreconact
                    pooldata.Sreconact = false(datasize, poolsize);
                    % or to use CSR (mostly in deployment)
                    pooldata.CSR = zeros(datasize, poolsize, 'int32');
                end
            case 'sharedC'
                if gpuDevice > 0
                    % CUDA
                    datapoolsize = uint64(datasize * databasis * poolsize * 4);
                    pooldata.image = calllib('CUpipeline', 'cudaworkspaceNew', datapoolsize); % single
                    SSRsize = uint64(datasize * poolsize);
                    pooldata.Sreconact =  calllib('CUpipeline', 'cudaworkspaceNew', SSRsize); % uint8
                    pooldata.CSR = calllib('CUpipeline', 'cudaworkspaceNew', SSRsize * 4); % int32
                else
                    % lib.pointer
                    pooldata.image = libpointer('singlePtr', zeros(datasize*databasis, poolsize));
                    pooldata.Sreconact =  libpointer('uint8Ptr', false(datasize, poolsize));
                    pooldata.Sreconact =  libpointer('int32Ptr', false(datasize, poolsize));
                end
            otherwise
        end
    otherwise
        1;
end
end
