function m = classsize(class)
% get Bytes of a class

% remove the Ptr
if endsWith(class, 'Ptr')
    class = class(1:end-3);
end

switch class
    case {'double', 'int64', 'uint64'}
        m = 8;
    case {'single', 'int32', 'uint32'}
        m = 4;
    case {'int24', 'uint24'}
        m = 3;
    case {'int16', 'uint16', 'half'}
        m = 2;
    case {'int8', 'uint8', 'char', 'logical', 'bool', 'string'}
        m = 1;
    case 'gpuArray'
        error('The gpuArray is illeagal in this function, plz replace the class() by classGPU() or gather the data to CPU.');
    otherwise
        m = 0;
end

end