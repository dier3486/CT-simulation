function varargout = gpuArraygroup(varargin)

for ii = 1:nargin
  varargout{ii} = gpuArray(varargin{ii});
end

end