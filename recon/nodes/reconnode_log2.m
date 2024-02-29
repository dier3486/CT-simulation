function [dataflow, prmflow, status] = reconnode_log2(dataflow, prmflow, status)
% recon node, log2
% [dataflow, prmflow, status] = reconnode_log2(dataflow, prmflow, status);

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

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);

% pipeline_onoff
if isfield(nodeprm, 'pipeline_onoff')
    pipeline_onoff = status.pipeline_onoff & nodeprm.pipeline_onoff;
else
    pipeline_onoff = status.pipeline_onoff;
end

% pipeline consol
dataflow_redirect = struct();
if pipeline_onoff
    statuscurr = status.pipepool.(nodename);
    nextnode = status.pipeline.(nodename).nextnode;
    % datasize
    datasize = max(0, statuscurr.WritePoint - statuscurr.ReadPoint);
    % writenum
    if ~isempty(nextnode)
        statusnext = status.pipepool.(nextnode);
        writenum = min(datasize, min(statusnext.WriteEnd, statusnext.poolsize) - statusnext.WritePoint + 1);
    else
        writenum = datasize;
        % Even when the nextnode is not existing the node will do its work as usual.
    end
    % copy pool to dataflow_redirect
    [dataflow_redirect, ~] = pooldatacopy(dataflow.pipepool.(nodename), ...
        dataflow_redirect, statuscurr.ReadPoint, 1, writenum, {}, 1);
else
    dataflow_redirect.rawdata = dataflow.rawdata;
    dataflow_redirect.rawhead = dataflow.rawhead;
end

% offset
if isfield(dataflow, 'offset')
    dataflow_redirect.offset = dataflow.offset;
    % Note: the offset will not be put in dataflow.pipepool.(nodename)!
end

% what it was
[dataflow_redirect, prmflow] = log2Kernelfuntion(dataflow_redirect, prmflow);

% copy back
if pipeline_onoff
    if ~isempty(nextnode)
        % copy dataflow_redirect to next pool
        [dataflow.pipepool.(nextnode), ~] = pooldatacopy(dataflow_redirect, ...
            dataflow.pipepool.(nextnode), 1, statusnext.WritePoint, writenum);
        % move next pool's write point
        status.pipepool.(nextnode).WritePoint = status.pipepool.(nextnode).WritePoint + writenum;
    end
    % move current pool's read point
    status.pipepool.(nodename).ReadPoint = status.pipepool.(nodename).ReadPoint + writenum;
    if datasize == 0
        % donothing did nothing, pass
        status.jobdone = 3;
    elseif writenum < datasize
        % done and keep waking
        status.jobdone = 2;
    else
        % normally done
        status.jobdone = 1;
    end
else
    dataflow.rawdata = dataflow_redirect.rawdata;
    status.jobdone = true;
end

status.errorcode = 0;
status.errormsg = [];
end


function [dataflow, prmflow] = log2Kernelfuntion(dataflow, prmflow)

% Z0
if isfield(prmflow, 'system') && isfield(prmflow.system, 'DBBzero')
    Z0 = prmflow.system.DBBzero;
else
    Z0 = 16384;
end
% offset
if isfield(dataflow, 'offset')
    dataflow.rawdata = dataflow.rawdata - mean(dataflow.offset.rawdata, 2);
else
    dataflow.rawdata = dataflow.rawdata - Z0;
end

Sneg = dataflow.rawdata<=0;
% if any(Sneg(:))
rawneg = -dataflow.rawdata(Sneg)./log(2) + 1/log(2);

% log2
dataflow.rawdata = -log2(dataflow.rawdata);
dataflow.rawdata(Sneg) = rawneg;

% intigration time
dataflow.rawdata = dataflow.rawdata + log2(single(dataflow.rawhead.Integration_Time));

end