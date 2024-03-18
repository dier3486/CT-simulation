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
% plstruct0.seriesdone = false;   % if true, this node will never be called in this series
% the status.pipeline.(nodes) is to control the states of the pipeline nodes; the consol node will (only) use these structs to
% control the pipeline work flow.

% loadrawdata is a default first node, need not to set in pipeline
if isfield(prmflow.protocol, 'loadrawdata')
    loadrawdata_onoff = prmflow.protocol.loadrawdata;
else
    loadrawdata_onoff = true;
end
% but user can close it by setting protocol.loadrawdata=false;

% initial status.pipeline and pipepool
status.pipeline = struct();
status.pipepool = struct();
% % we will delete the status.pipepool, later

% pine-line node 'loadrawdata'
if loadrawdata_onoff
    status.pipeline.loadrawdata = plstruct0;
    status.pipeline.loadrawdata.sleeping = false;   % wake up the node loadrawdata
    status.pipepool.loadrawdata = struct(); 
end

% initial other pine-line nodes
pipenodes = fieldnames(prmflow.pipe);
if ~isempty(pipenodes)
    % nextnode of loadrawdata
    if loadrawdata_onoff
        status.pipeline.loadrawdata.nextnode = pipenodes{1};
    end
    % other nodes
    for ii = 1:length(pipenodes)
        % ini status.pipeline.(pipenodes{ii})
        status.pipeline.(pipenodes{ii}) = plstruct0;
        status.pipeline.(pipenodes{ii}).index = ii;
        status.pipeline.(pipenodes{ii}).priority = ii;
        if ii < length(pipenodes)
            status.pipeline.(pipenodes{ii}).nextnode = pipenodes{ii + 1};
        end
        if ii > 1
            status.pipeline.(pipenodes{ii}).prevnode = pipenodes{ii - 1};
        else
            if loadrawdata_onoff
                status.pipeline.(pipenodes{ii}).prevnode = 'loadrawdata';
            else
                status.pipeline.(pipenodes{ii}).prevnode = 'NULL';
            end
        end
        % pipeline_onoff
        if isfield(prmflow.pipe.(pipenodes{ii}), 'pipeline_onoff')
            status.pipeline.(pipenodes{ii}).pipeline_onoff = ...
                status.pipeline.(pipenodes{ii}).pipeline_onoff & prmflow.pipe.(pipenodes{ii}).pipeline_onoff;
        end
        % ini status.pipepool
        status.pipepool.(pipenodes{ii}) = struct();
    end
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
status.defaultpooldata = struct();
status.defaultpooldata.rawhead = struct();
status.defaultpooldata.rawdata = single([]);
% % We shall make it configuerable later.

% default public buffer fields
status.defaultpublicpool = setdefaultpublicpool(status.defaultpooldata);

% we will move these to dataflow
pipepoolnodes = fieldnames(status.pipepool);
for ii = 1:length(pipepoolnodes)
    status.pipepool.(pipepoolnodes{ii}) = status.defaultpublicpool;
end
% Those default settings can be changed in node's prepare.

% initial pipepool and buffer in dataflow
dataflow.pipepool = struct();
dataflow.buffer = struct();

% NULL pool
status.pipepool.NULL = [];
dataflow.pipepool.NULL = [];

% default private buffer fields for the nodes (suggest)
status.defaultprivatepool = setdefaultprivatepool(status.defaultpooldata);
% to use it in nodes' prepare like this: 
%   dataflow.buffer.(nodename).outpool = status.defaultprivatepool;
% Anyhow, it is not forced, the nodes' owners are permitted to use any structure format in whose private buffer.

end


function publicpool = setdefaultpublicpool(defaultpooldata)

% default public buffer fields
publicpool = struct();
% the point for writing data
publicpool.WritePoint = 1;
% the end point in planing of the written data (use to limit the written data size)
publicpool.WriteEnd = Inf;
% the point for reading data
publicpool.ReadPoint = 1;
publicpool.ReadEnd = Inf;
publicpool.AvailNumber = 0;     % AvailNumber==WritePoint-ReadPoint.
% poolsize
publicpool.poolsize = Inf;
% Normally the written data range is [WritePoint, min(WritePoint-1, min(WriteEnd, poolsize)].
publicpool.warningstage = 1024;
publicpool.recylestrategy = 1;
% 0: only when filled, 1: always, 2: over warningstage.

publicpool.datafields = fieldnames(defaultpooldata);
publicpool.data = defaultpooldata;

end

function privatepool = setdefaultprivatepool(defaultpooldata)
% default private buffer fields for the nodes (suggest)
privatepool = struct();

privatepool.ReadPoint = 1;
privatepool.WritePoint = 1;
privatepool.WriteEnd = Inf;
privatepool.ReadEnd = Inf;
privatepool.AvailNumber = 0;
privatepool.poolsize = Inf;
privatepool.recylestrategy = 1;
% 0: only when filled (for axial), 1: always.
% extra
privatepool.circulatemode = false;      % circulation buffer
privatepool.ReadViewindex = 1;          % Viewindex if the ReadPoint
% privatepool.AvailViewindex = 0;         % Viewindex if the AvailPoint
% privatepool.ReadStuck = false;          % onoff to close the reading
% privatepool.WriteStuck = false;         % onoff to close the writing

privatepool.datafields = fieldnames(defaultpooldata);
privatepool.data = defaultpooldata;

end