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

% pine-line node 'loadrawdata'
status.pipeline = struct();
status.pipeline.loadrawdata = plstruct0;
status.pipeline.loadrawdata.sleeping = false;   % wake up the node loadrawdata
% % the prepare of loadrawdata
% if isempty(prmflow.protocol.datablock)
%     % no blocks
%     status.buffer.loadrawdata.datablock_onoff = false;
%     status.buffer.loadrawdata.Nblock = 1;
%     status.buffer.loadrawdata.datablocksize = prmflow.protocol.viewnumber;
% else
%     status.buffer.loadrawdata.datablock_onoff = true;
%     viewnumber = prmflow.protocol.viewnumber;
%     datablock = prmflow.protocol.datablock;
%     Nblock = ceil(viewnumber / datablock);
%     status.buffer.loadrawdata.Nblock = Nblock;
%     status.buffer.loadrawdata.datablocksize = ...
%         [repmat(datablock, 1, Nblock-1)  viewnumber - datablock * (Nblock-1)];
% end
% status.buffer.loadrawdata.iblock = 1;
% % to prepare a 'hidden' node of loadrawdata in pipe-line
% status.buffer.loadrawdata.pipepool = struct();
% status.buffer.loadrawdata.pipepool.ReadPoint = 1;
% status.buffer.loadrawdata.pipepool.WritePoint = 1;
% % hard-codes (we need a readrawdataprepare)
% status.buffer.loadrawdata.pipepool.poolsize = Inf;
% status.buffer.loadrawdata.pipepool.warningstage = 2048;
% status.buffer.loadrawdata.pipepool.recylestrategy = 0;
% % The 'loadrawdata' is a forced node, the prepare of loadrawdata is also forced. If user configured another node to read
% % rawdata (by calling reconnode_readrawdata), the prepare function reconnode_readrawdataprepare (if exist) will be called to
% % do the prepare for that user configured node. They are different things.

% pipepool 'loadrawdata'
status.pipepool = struct();
status.pipepool.loadrawdata = struct(); 

% initial status.pipeline and pipepool
pipenodes = fieldnames(prmflow.pipe);
if ~isempty(pipenodes)
    % nextnode of loadrawdata
    status.pipeline.loadrawdata.nextnode = pipenodes{1};
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
            status.pipeline.(pipenodes{ii}).prevnode = 'loadrawdata';
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
    % Even when the pipeline_onoff is false, there still have a status.pipeline which includes all the pineline nodes.
end

% default pool fields are rawhead, rawdata and to be set up
status.defaultpool = struct();
status.defaultpool.rawhead = struct();
status.defaultpool.rawdata = single([]);
% % We shall make it configuerable later.


% I know the status.pipepool is the I/O buffer of the pipeline nodes; 
% the status.buffer is the private buffers of the pipeline nodes.

% orgnodename = status.nodename;
% % loop the pipenodes
% if ~isempty(pipenodes)
%     for ii = 1:length(pipenodes)
%         % ini status.buffer.(pipenodes{ii}) and status.pipepool.(pipenodes{ii}) = struct();
% %         status.buffer.(pipenodes{ii}) = struct();
%         status.pipepool.(pipenodes{ii}) = struct();
%         % nodes initial function
%         nodename_slip = regexp(pipenodes{ii}, '_', 'split');
%         inifunname = ['reconnode_' lower(nodename_slip{1}) 'initial'];
%         % but we have not any nodes initial functions yet 
%         status.nodename = pipenodes{ii};
%         if any(exist(inifunname) == [2 5 6])
%             inifun = str2func(inifunname);
%             [prmflow, status] = inifun(prmflow, status);
%             % the status.buffer and status.pipepool should be done in the prepare function.
%         end
% %         else
% %             % default prepare of pipepool
% %             status.pipepool.(pipenodes{ii}).rawhead = struct();
% %             status.pipepool.(pipenodes{ii}).rawdata = single([]);
% %             for jj = 1:length(status.defaultpool)
% %                 status.pipepool.(pipenodes{ii}).(status.defaultpool{jj}) = cast([], 'single');
% %             end
% %         end
%     end
% end
% status.nodename = orgnodename;

% The 'preparefun's have put the fields in status.pipepool and/or status.buffer, which they want to occur in the
% dataflow.pipepool and/or dataflow.buffer.
% DO NOT put data in status.pipepool or status.buffer, they are configures and/or initial states.
% In runing the pipeline the data shall be put in dataflow.

% we don't like these variables exist in dataflow
pipepoolnodes = fieldnames(status.pipepool);
for ii = 1:length(pipepoolnodes)
    % the point for writing data
    status.pipepool.(pipepoolnodes{ii}).WritePoint = 1;
    % the end point in planing of the written data (use to limit the written data size)
    status.pipepool.(pipepoolnodes{ii}).WriteEnd = Inf;
    % Normally the written data range is [WritePoint, min(WritePoint-1, WriteEnd)].
    % the point for reading data
    status.pipepool.(pipepoolnodes{ii}).ReadPoint = 1;
    status.pipepool.(pipepoolnodes{ii}).ReadEnd = Inf;
    % We have not used the WriteEnd and ReadEnd yet.
    if ~isfield(status.pipepool.(pipepoolnodes{ii}), 'poolsize')
        status.pipepool.(pipepoolnodes{ii}).poolsize = Inf;
    end
    if ~isfield(status.pipepool.(pipepoolnodes{ii}), 'warningstage')
        status.pipepool.(pipepoolnodes{ii}).warningstage = 1024;
    end
    if ~isfield(status.pipepool.(pipepoolnodes{ii}), 'recylestrategy')
        status.pipepool.(pipepoolnodes{ii}).recylestrategy = 1;
        % 0: never, 1: always, 2: over warningstage.
    end
    % flag to sign if all of the raw data are ready, that not new data will be written in. (not used)  
    status.pipepool.(pipepoolnodes{ii}).EndData = false;
end

% initial pipepool and buffer in dataflow
dataflow.pipepool = struct();
dataflow.buffer = struct();

end