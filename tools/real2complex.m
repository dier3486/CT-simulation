function varargout = real2complex(varargin)
% if x is real, x = x+i*x; else x = x; and so on
% [x, y, z] = real2complex(x, y, z);

s = cellfun(@isreal, varargin);
varargout = varargin;
varargout(s) = cellfun(@complex, varargout(s), varargout(s), 'UniformOutput', false);

end