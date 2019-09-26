function cfgstruct = materialconfigure(cfgstruct, samplekeV, elementspath)
% configure the materials in a struct

if nargin < 2
    samplekeV = [];
end
if nargin < 3
    elementspath = '';
end

Ncfg = length(cfgstruct(:));
if iscell(cfgstruct)
    % recurse
    for icfg = 1:Ncfg
        cfgstruct{icfg} = materialconfigure(cfgstruct{icfg}, samplekeV, elementspath);
    end
elseif isstruct(cfgstruct)
    if Ncfg>1
        % recurse
        for icfg = 1:Ncfg
            cfgstruct(icfg) = materialconfigure(cfgstruct(icfg), samplekeV, elementspath);
        end
    else
        fields = fieldnames(cfgstruct);
        % recurse
        for ifield = 1:length(fields)
            if iscell(cfgstruct.(fields{ifield})) || isstruct(cfgstruct.(fields{ifield}))
                cfgstruct.(fields{ifield}) = ...
                    materialconfigure(cfgstruct.(fields{ifield}), samplekeV, elementspath);
            end
        end
        % define material
        if isfield(cfgstruct, 'material')
            if ischar(cfgstruct.material)
                % read configure file
                cfgstruct.material = loadmaterial(cfgstruct.material);
            end
            if isstruct(cfgstruct.material)
                cfgstruct.material = ...
                    materialdefine(cfgstruct.material, samplekeV, elementspath);
            end
        end
    end
end


