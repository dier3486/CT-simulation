function [dataflow, prmflow, status] = reconnode_FanRadialrebin(dataflow, prmflow, status)
% recon node, fan-radial rebin, the radial rebin before Azi
% [dataflow, prmflow, status] = reconnode_FanRadialrebin(dataflow, prmflow, status);
% Support (X)DFS

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

% nodename
nodename = status.nodename;

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% pipeline consol
dataflow_redirect = struct();
if pipeline_onoff
    statuscurr = status.pipepool.(nodename);
    nextnode = status.pipeline.(nodename).nextnode;
    % datasize
    datasize = max(0, statuscurr.WritePoint - statuscurr.ReadPoint);
    % only even
    datasize = floor(datasize/2)*2;
    % writenum
    if ~isempty(nextnode)
        statusnext = status.pipepool.(nextnode);
        writenum = min(datasize, min(statusnext.WriteEnd, statusnext.poolsize) - statusnext.WritePoint + 1);
        % only even
        writenum = floor(writenum/2)*2;
    else
        writenum = datasize;
        % Even when the nextnode is not existing the node will do its work as usual.
    end
    % copy pool to dataflow_redirect
    [dataflow_redirect, ~] = pooldatacopy(dataflow.pipepool.(nodename), ...
        dataflow_redirect, statuscurr.ReadPoint, 1, writenum, {}, 1);
end

% what it was
if pipeline_onoff
    [dataflow_redirect, prmflow] = FanRadialrebinKernelfunction(dataflow_redirect, prmflow, status);
else
    [dataflow, prmflow] = FanRadialrebinKernelfunction(dataflow, prmflow, status);
end

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
    status.jobdone = true;
end


% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


