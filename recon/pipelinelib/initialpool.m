function currpool = initialpool(currpool, nodepipeprm, namekey)
% initial a pool by user setting
% While the namekey is 'public' the currpool is supposed the dataflow.pipepool.(nodename), which will be used in
% defaultpipelineprepare. Or, any other namekey can be employed in the node's preapre.

poolfields = [namekey 'fields'];
poolclasses = [namekey 'classes'];
poolpoints = [namekey 'pool'];

if isfield(nodepipeprm, poolfields)
    poolfields = nodepipeprm.(poolfields);
    if ~iscell(poolfields)
        poolfields = splitcomma(poolfields);
    end
    if isfield(nodepipeprm, poolclasses)
        poolclass = nodepipeprm.(poolclasses);
        if ~iscell(poolclass)
            poolclass = splitcomma(poolclass);
        end
    else
        poolclass = {};
    end
else
    poolfields = {};
end

% ini data
if ~isempty(poolfields)
    currdata = struct();
    for ii = 1:length(poolfields)
        if ~isempty(poolclass)
            currdata.(poolfields{ii}) = uint8cast([], poolclass{ii});
        else
            currdata.(poolfields{ii}) = single([]);
        end
    end
    currpool.datafields = poolfields;
    currpool.data = currdata;
end

% poolpoints user configured
if isfield(nodepipeprm, poolpoints)
    poolpoints = nodepipeprm.(poolpoints);
else
    poolpoints = struct();
end
currpool = structmerge(poolpoints, currpool);
% Note: the user's configure will overwrite all above settings.

end