function y = castuint8(x, type)
% cast variable to uint8
%   y = castuint8(x);
% or, e.g.
%   y = castuint8(x, 'uint24');


% the 'type' is a modifier of the forced type-cast of x to 'type' before casting to uint8
if nargin < 2
    type = '';
    % empty means to skip it, which is suggested while the class(x)==type.
    % Due to not all the classes are supported in matlab, e.g. uint24, that
    % we need a modifier 'type' to present those extra classed.
end
type = lower(type);

% no GPU
if isa(x, 'gpuArray')
    x = gather(x);
end

% to keep the shape
sizeX = num2cell(size(x));
x = x(:);

% type-cast to uint8
if isnumeric(x)
    if ~isreal(x)
        x = reshape([real(x) imag(x)]', [], 1);
    end
    switch type
        case ''
            % default
            y = typecast(x, 'uint8');
        case {'float', 'float32', 'single'}
            % float is single
            y = typecast(cast(x, 'single'), 'uint8');
        case {'float64', 'double'}
            % double
            y = typecast(cast(x, 'double'), 'uint8');
        case 'uint24'
            % special type for 24bit rawdata
            y = reshape(typecast(cast(x, 'uint32'), 'uint8'), 4, []);
            s_overflow = y(4,:)>0;
            y = y(1:3, :);
            y(:, s_overflow) = uint8(255);
            y = y(:);
        case {'int64b', 'uint64b', 'int32b', 'uint32b', 'int16b', 'uint16b', 'int8b', 'uint8b', 'uint24b'}
            % big endian int
            y = castuint8(x, type(1:end-1));
            m = classsize(type(1:end-1));
            y = flipud(reshape(y, m, []));
            y = y(:);
        case 'uint8'
            % just cast
            y = cast(x, 'uint8');
        otherwise
            % other? 
            y = typecast(cast(x, type), 'uint8');
    end
else 
    switch class(x)
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
            % unkown class
            y = [];
    end
end

y = reshape(y, [], sizeX{2:end});

end