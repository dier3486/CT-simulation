function y = castuint8(x, type)
% cast viarible to uint8

if nargin < 2
    type = '';
end
switch class(x)
    case {'double', 'single', 'int64', 'uint64', 'int32', 'uint32', ...
          'int16', 'uint16', 'int8', 'uint8'}
        if isempty(type)
            y = typecast(x, 'uint8');
        else
            y = typecast(cast(x, type), 'uint8');
        end
    case 'char'
        % y = typecast(uint16(x), 'uint8');
        y = uint8(x);
    case 'logical'
        % y = uint8(bin2dec(num2str(rot90(reshape(x, 8, []),-1))));
        y = uint8(x);
    case 'cell'
        tmp = cell(size(x));
        for ii = 1:length(x)
            tmp{ii} = castuint8(x{ii}, type);
        end
        y = cell2mat(tmp);
    otherwise
        y = [];
end