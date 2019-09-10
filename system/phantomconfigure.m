function phantom = phantomconfigure(phantom_cfg)
% configure a phantom

if ~isstruct(phantom_cfg)
    % read configure file
    if exist(phantom_cfg, 'file')
        [~, ~, cfgext] = fileparts(phantom_cfg);
        
    elseif 1
    end
        
    
    phantom_cfg = [];
end

% copy
phantom = phantom_cfg;
% object to cell
if ~iscell(phantom.object)
    phantom.object = num2cell(phantom.object);
end

