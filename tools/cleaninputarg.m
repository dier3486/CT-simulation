function varargout = cleaninputarg(defaultval, varargin)
% clean input arguments for functions
% assume a function myfun(x, y, varargin) have default values of the varargin {a0, b0, c0}, call this
% [a, b, c] = cleaninputarg({a0, b0, c0}, varargin{:});
% to fill up the empty and skipped input arguments
% NOTE: here the input is the 'varargin{:}' but NOT the 'varargin'! Nor it will be confused in a cell argument or a 
% multi-argument input. 

if ~iscell(defaultval)
    defaultval = {defaultval};
end

N = min(length(defaultval), nargin-1);
s = ~cellfun(@isempty, varargin(1:N));
defaultval(s) = varargin(s);
varargout = defaultval(1:nargout);

end