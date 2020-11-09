function [GPUinfo, errinfo] = initialGPU(index)
% initialGPU

if nargin<1
    index = 1;
end

errinfo = [];
if index
    try
        GPUinfo = gpuDevice;
    catch me
        % no GPU
        GPUinfo = [];
        errinfo = me;
        return
    end
    if GPUinfo.Index ~= index
        % reselect GPU device
        GPUinfo = gpuDevice(index);
    end
else
    GPUinfo = [];
end

end