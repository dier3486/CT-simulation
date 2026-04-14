function varargout = putinGPU(varargin)
% put the inputs in GPU, like
%  [a, b, c] = putinGPU(a, b, c);

for ii = 1:nargin
    varargout{ii} = everything2single(varargin{ii}, 'any', 'gpuArray');
end

% NOTE: to gather things from GPU inversely call mystruct = everything2single(mystruct, 'any', 'gather');
end