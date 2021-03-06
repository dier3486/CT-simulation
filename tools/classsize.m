function m = classsize(class)
% get Bytes of a class

switch class
    case {'double', 'int64', 'uint64'}
        m = 8;
    case {'single', 'int32', 'uint32'}
        m = 4;
    case {'int24', 'uint24'}
        m = 3;
    case {'int16', 'uint16'}
        m = 2;
    case {'int8', 'uint8', 'char', 'logical', 'bool'}
        m = 1;
    otherwise
        m = 0;
end

end