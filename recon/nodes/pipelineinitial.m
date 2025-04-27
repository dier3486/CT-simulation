function [dataflow, prmflow, status] = pipelineinitial(dataflow, prmflow, status)
% subfunction in recon initial for data-blocks and pipe-line
% [dataflow, prmflow, status] = pipelineinitial(dataflow, prmflow, status)

% Copyright Dier Zhang
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% pipeline_onoff
status.pipeline_onoff = prmflow.protocol.pipelinereplicate;

% ini pipeline.(nodes)
plstruct0 = struct();
plstruct0.index = 0;            % node index
plstruct0.nextnode = '';        % next node
plstruct0.prevnode = '';        % previous node
plstruct0.priority = 0;         % node priority
plstruct0.sleeping = true;      % node running status
plstruct0.prepared = false;     % node prepared flag
plstruct0.pipeline_onoff = status.pipeline_onoff;   % node pipeline_onoff flag
plstruct0.walltime = 0;
plstruct0.jobdone = 0;          % to record the previous jobdone flag
plstruct0.branchnextnodes = {};         % branch nodes
plstruct0.branchnextpoolindex = [];
% the status.pipeline.(nodes) is to control the states of the pipeline nodes; the consol node will (only) use these structs to
% control the pipeline work flow.


% loadrawdata is a default first node, need not to set in pipeline
if isfield(prmflow.protocol, 'loadrawdata')
    loadrawdata_onoff = prmflow.protocol.loadrawdata;
else
    loadrawdata_onoff = true;
end
% but user can close it by setting protocol.loadrawdata=false;

% initial status.pipeline
status.pipeline = struct();

% pine-line node 'loadrawdata'
if loadrawdata_onoff && ~isfield(prmflow.pipe, 'loadrawdata')
    % add the loadrawdata to prmflow.pipe
    prmflow.pipe.loadrawdata = struct();
    % move it to first
    Nnode = length(fieldnames(prmflow.pipe));
    prmflow.pipe = orderfields(prmflow.pipe, [Nnode 1:Nnode-1]);
end

% iteration or other looping nodes
prmflow.pipe =  combinenodes(prmflow.pipe);

% initial other pine-line nodes
pipenodes = fieldnames(prmflow.pipe);
if ~isempty(pipenodes)
    for ii = 1:length(pipenodes)
        % ini status.pipeline.(pipenodes{ii})
        status.pipeline.(pipenodes{ii}) = plstruct0;
        status.pipeline.(pipenodes{ii}).index = ii;
        status.pipeline.(pipenodes{ii}).priority = ii;
        % next node
        if ii < length(pipenodes)
            status.pipeline.(pipenodes{ii}).nextnode = pipenodes{ii + 1};
        else
            status.pipeline.(pipenodes{ii}).nextnode = 'NULL';
        end
        % previous node
        if ii > 1
            status.pipeline.(pipenodes{ii}).prevnode = pipenodes{ii - 1};
        else
            status.pipeline.(pipenodes{ii}).prevnode = 'NULL';
        end
        % pipeline_onoff
        if isfield(prmflow.pipe.(pipenodes{ii}), 'pipeline_onoff')
            status.pipeline.(pipenodes{ii}).pipeline_onoff = ...
                status.pipeline.(pipenodes{ii}).pipeline_onoff & prmflow.pipe.(pipenodes{ii}).pipeline_onoff;
        end
    end
    % to wake up the first node
    status.pipeline.(pipenodes{1}).sleeping = false;
end

% if ~pipeline_onoff, return
if ~prmflow.protocol.pipelinereplicate
    % pipeline off, but to read data in blocks.
    if ~isempty(prmflow.protocol.datablock)
        status.pipeline.loadrawdata.priority = 999;
        % It is a debug mod, loadrawdata will have highest priority to loop the blocks until all data being read out before the
        % next node waking up, which equivlaent to non-block (1 block) data reading.
    end
    % job done whether ~pipeline_onoff
    return;
    % Note : even when the pipeline_onoff is false, there still have a status.pipeline which includes all the pineline nodes.
end

% default pool fields are rawhead, rawdata and to be set up
defaultpooldata = struct();
defaultpooldata.rawhead = struct();
defaultpooldata.rawdata = single([]);
% We shall make it configuerable later.

% default public buffer fields
status.defaultpool = setdefaultpool(defaultpooldata);

% initial pipepool and buffer in dataflow
dataflow.pipepool = struct();
dataflow.buffer = struct();

% NULL pool
dataflow.pipepool.NULL = [];

% % default private buffer fields for the nodes (suggest)
% status.defaultprivatepool = setdefaultpool(status.defaultpooldata);
% to use it in nodes' prepare like this: 
%   dataflow.buffer.(nodename).outpool = status.defaultprivatepool;
% Anyhow, it is not forced, the nodes' owners are permitted to use any structure format in whose private buffer.

end


function defaultpool = setdefaultpool(defaultpooldata)
% default private buffer fields for the nodes (suggest)
defaultpool = struct();

% poolsize
defaultpool.poolsize = Inf;
% the point for writing data
defaultpool.WritePoint = 1;

% the start/end point in writing data, which is the first/last writable position.
defaultpool.WriteStart = 1;
defaultpool.WriteEnd = Inf;
% Normally the wrottable buffer range is [WritePoint, min(WriteEnd, poolsize)].

% the point for reading data
defaultpool.ReadPoint = 1;

% the start/end point in reading, which is the first/last readable position. 
defaultpool.ReadStart = 1;
defaultpool.ReadEnd = Inf;
% the ReadStart and ReadEnd are mostly used to asign the start and end of a shot as a shot lock 

% the available point
defaultpool.AvailPoint = 0;

% % the available data size
% defaultpool.AvailNumber = 0;    

% recyling strategy
defaultpool.recylestrategy = 1;
% 0: never (*), 1: always, 2: over warningstage.
% *: but for a circulated buffer it could be clear to 0 in necessary moments.
% the warningstage
defaultpool.warningstage = Inf;
% to keep some bottom in recyling
defaultpool.keepbottom = 0;

% prevWritePoint
defaultpool.prevWritePoint = 1;
% use to record the WritePoint before the writting

% circulated buffer, onoff
defaultpool.circulatemode = false;
% we can employ a circulated buffer for axial to simplify the periodic algorithm.

% flag isshotstart
defaultpool.isshotstart = true;
% use to judge if to initialize the pool, better than to use WritePoint==WriteStart.

% % shot end point
% defaultpool.ShotEnd = Inf;
% the function of the ShotEnd is covered by/ the ReadEnd.

% Viewindex if the ReadPoint
defaultpool.ReadViewindex = 1;

% Viewindex if the AvailPoint (ReadPoint + AvailNumber - 1)
defaultpool.AvailViewindex = 0;     

% onoff to close the writing and/or reading
defaultpool.WriteStuck = false;
defaultpool.ReadStuck = false;

% the buffer of data in the pool
defaultpool.datafields = fieldnames(defaultpooldata);
defaultpool.data = defaultpooldata;
% Note: Some times it is just a declaration, the real buffer could be in other space, especially when we put the pool in root
% structure 'status' which is not suitable to handle a buffer resource.

% is carried and carry pool
defaultpool.iscarried = false;
defaultpool.carrynode = '';
defaultpool.carryindex = 1;  % the carrynode's inputpool's index
% if carried the defaultpool.data shall be a null struct that the real data buffer is not in the current pool,
% which shall be found in dataflow.pipepool.(carrynode)(carryindex).data.
defaultpool.carriages = {};

% buffer resource
defaultpool.bufferresource = '';
% It is a label of the buffer resource, e.g. 'GPU device1', not used yet.

% trace record
defaultpool.trace(1) = poolmirror(defaultpool);
defaultpool.trace(1).operator = 'default';
% in debug mode, the trace can be used to record those points and flags as a log.



end