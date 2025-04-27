function hout = cttool(varargin)
% entry

hFig = imcttool(varargin{:});
if (nargout > 0)
    % Only return handle if caller requested it.
    hout = hFig;
end

end
