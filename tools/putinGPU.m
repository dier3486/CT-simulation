function varargout = putinGPU(varargin)
% put the inputs in GPU

for ii = 1:nargin
    varargout{ii} = everything2single(varargin{ii}, 'any', 'gpuArray');
end

end