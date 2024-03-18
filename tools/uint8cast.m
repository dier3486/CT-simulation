function y = uint8cast(x, type)
% cast uint8 to variable 

switch lower(type)
    case {'double', 'single', 'int64', 'uint64', 'int32', 'uint32', ...
          'int16', 'uint16', 'int8', 'uint8'}
        y = typecast(x, type);
    case {'float', 'float32'}
        y = typecast(x, 'single');
    case 'char'
        % y = char(typecast(x, 'uint16'));
        y = char(x);
    case {'logical', 'bool'}
        % y = rot90(dec2bin(x, 8)=='1');
        y = logical(x);
    case {'uint24'}
        % special cast type for 24bit rawdata
        x = reshape(x, 3, []);
        x = [x; zeros(1, size(x, 2))];
        y = typecast(x(:), 'uint32');
    case {'int64b', 'uint64b', 'int32b', 'uint32b', 'int16b', 'uint16b', 'int8b', 'uint8b', 'uint24b'}
        % big endian int
        m = classsize(type(1:end-1));
        x = flipud(reshape(x, m, []));       
        y = uint8cast(x(:), type(1:end-1));
        
    % other special cases (may not used in sparsepack)
    case 'struct'
        % an empty struct
        y = struct();
    case 'cell'
        % to cell
        y = num2cell(x);
    otherwise
        % empty
        y = [];
end

end