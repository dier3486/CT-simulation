function r = classGPU(x)
% I can not remember what classUnderlying nor underlyingtype

switch class(x)
    case 'gpuArray'
        r = classUnderlying(x);
    otherwise
        r = class(x);
end