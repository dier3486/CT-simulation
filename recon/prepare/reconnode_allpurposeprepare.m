function [dataflow, prmflow, status] = reconnode_allpurposeprepare(dataflow, prmflow, status)
% prepare node, prepare to do nothing
% [dataflow, prmflow, status] = reconnode_allpurposeprepare(dataflow, prmflow, status);
% wow~

% parameters set in pipe
nodename = status.nodename;
if isfield(prmflow.pipe, nodename)
    nodeprm = prmflow.pipe.(nodename);
else
    nodeprm = struct();
end

% node level
if isfield(nodeprm, 'level')
    nodelevel = nodeprm.level;
else
    nodelevel = 0;
    prmflow.pipe.(nodename).level = 0;
end

% main operator
if isfield(nodeprm, 'mainopt')
    if ischar(nodeprm.mainopt)
        prmflow.pipe.(nodename).mainopt = eval(nodeprm.mainopt);
    end
else
    prmflow.pipe.(nodename).mainopt = @(x) x+1;
end

% echo onoff
if ~isfield(nodeprm, 'pipeinfoecho')
    prmflow.pipe.(nodename).pipeinfoecho = true;
end


% node function
if isfield(nodeprm, 'funct')
    nodefunct = nodeprm.funct;
else
    prmflow.pipe.(nodename).funct = '';
    nodefunct = '';
end

% prepare by existing nodes
if ~isempty(nodefunct)
    preparename = ['reconnode_' lower(nodefunct) 'prepare'];
    if any(exist(preparename) == [2 5 6]) % can run
        if status.echo_onoff
            fprintf(['(' nodefunct ')']);
        end
        preparefun = str2func(preparename);
        [dataflow, prmflow, status] = preparefun(dataflow, prmflow, status);
        return;
    end
end


% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

if pipeline_onoff
    % public input, (nodeprm.inputfields)
    [status.pipepool.(nodename), dataflow.pipepool.(nodename)] = ...
        initialpool(status.pipepool.(nodename), nodeprm, 'public', 'public', status);

    % private output, (nodeprm.outputfields)
    if nodelevel > 0
        [dataflow.buffer.(nodename).outputpool, ~] = ...
            initialpool(status.pipepool.(nodename), nodeprm, 'output', 'private', status);
    end

    % private input
    if nodelevel > 1
        [dataflow.buffer.(nodename).inputpool, ~] = ...
            initialpool(status.pipepool.(nodename), nodeprm, 'input', 'private', status);
    end
    % private inner
    if nodelevel > 2
        [dataflow.buffer.(nodename).innerpool, ~] = ...
            initialpool(status.pipepool.(nodename), nodeprm, 'inner', 'private', status);
    end
end

status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end


function [currpool, currdata] = initialpool(currpool, nodeprm, namekey, typeflag, status)
% currpool is status.pipepool.(nodename), currdata is dataflow.pipepool.(nodename) when typeflag=='public';
% Or, currpool is dataflow.buffer.(nodename).([namekey 'pool']), currdata is ~ when typeflag=='private'.

poolfields = [namekey 'fields'];
poolclasses = [namekey 'classes'];
poolpoints = [namekey 'pool'];

if isfield(nodeprm, poolfields)
    poolfields = nodeprm.(poolfields);
    if ~iscell(poolfields)
        poolfields = splitcomma(poolfields);
    end
    if isfield(nodeprm, poolclasses)
        poolclass = nodeprm.(poolclasses);
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
if isempty(poolfields)
    currdata = status.defaultpooldata;
else
    currdata = struct();
    for ii = 1:length(poolfields)
        if ~isempty(poolclass)
            currdata.(poolfields{ii}) = uint8cast([], poolclass{ii});
        else
            currdata.(poolfields{ii}) = single([]);
        end
    end
end

% ini pool
switch typeflag
    case {'public', 'publ'}
        currpool = status.defaultpublicpool;
        currpool.datafields = poolfields;
        % I know pool data is not in the public pool (We should do that)
    case {'private', 'priv'}
        currpool = status.defaultprivatepool;
        currpool.datafields = poolfields;
        currpool.data = currdata;
    otherwise
        % error
        error(0);
end

% poolpoints
if isfield(nodeprm, poolpoints)
    poolpoints = nodeprm.(poolpoints);
else
    poolpoints = struct();
end
currpool = structmerge(poolpoints, currpool);

end