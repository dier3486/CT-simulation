function phantom = phantomconfigure(phantom_cfg)
% configure a phantom

% % debug
% phantom_cfg = './test phantom.xml';

if ~isstruct(phantom_cfg)
    if ischar(phantom_cfg)
        % read configure file
        phantom_cfg = readcfgfile(phantom_cfg);
    else
        error(['Illegal configure class ', class(phantom_cfg)]);
    end
end

% copy
phantom = phantom_cfg;
% object to cell
if ~iscell(phantom.object)
    phantom.object = num2cell(phantom.object);
end
% fill up default fields
phantom = defaultphantom(phantom);

end


function phantom = defaultphantom(phantom)
% default fields of a phantom

% Nobject
if ~isfield(phantom, 'Nobject')
    phantom.Nobject = length(phantom.object(:));
end
% object_root
if ~isfield(phantom, 'object_root')
    phantom.object_root = 0;
end
% object_tree
if ~isfield(phantom, 'object_tree')
    phantom.object_tree = zeros(1, phantom.Nobject);
end
% objects
for iobj = 1:phantom.Nobject
    % vector
    phantom.object{iobj}.vector = reshape(phantom.object{iobj}.vector, 3, 3);
    % invV
    if ~isfield(phantom.object{iobj}, 'invV')
        phantom.object{iobj}.invV = inv(phantom.object{iobj}.vector);
    end
    % volume
    if ~isfield(phantom.object{iobj}, 'volume')
        phantom.object{iobj}.volume = ...
            defaultvolume(phantom.object{iobj}.vector, phantom.object{iobj}.type);
    end
end

end
